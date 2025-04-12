use std::path::Path;
use std::io::Write;
use std::collections::HashSet;
use anyhow::Result;
use noodles::fastq::{self, Record};
use noodles::fastq::AsyncReader;
use twox_hash::XxHash64;
use bstr::ByteSlice;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader, AsyncWrite, AsyncWriteExt};
use tokio::sync::{mpsc, Mutex as TokioMutex};
use tokio::task;
use futures::{Stream, StreamExt};
use std::sync::{Arc, Mutex};

// Number of records to process in a batch
const BATCH_SIZE: usize = 5000;

#[derive(Debug, serde::Serialize, serde::Deserialize, Clone)]
pub struct FastqDeduplicationStats {
    total_reads: u64,
    unique_reads: u64,
    duplicate_reads: u64,
}

impl FastqDeduplicationStats {
    fn new() -> Self {
        Self {
            total_reads: 0,
            unique_reads: 0,
            duplicate_reads: 0,
        }
    }

    fn increment_total(&mut self) {
        self.total_reads += 1;
    }

    fn increment_unique(&mut self) {
        self.unique_reads += 1;
    }

    fn increment_duplicate(&mut self) {
        self.duplicate_reads += 1;
    }

    fn merge(&mut self, other: &Self) {
        self.total_reads += other.total_reads;
        self.unique_reads += other.unique_reads;
        self.duplicate_reads += other.duplicate_reads;
    }

    fn print(&self) {
        println!("Total reads: {}", self.total_reads);
        println!("Unique reads: {}", self.unique_reads);
        println!("Duplicate reads: {}", self.duplicate_reads);
    }
}

impl IntoPy<PyObject> for FastqDeduplicationStats {
    fn into_py(self, py: Python) -> PyObject {
        let dict = PyDict::new_bound(py);
        dict.set_item("total_reads", self.total_reads).unwrap();
        dict.set_item("unique_reads", self.unique_reads).unwrap();
        dict.set_item("duplicate_reads", self.duplicate_reads).unwrap();
        dict.into()
    }
}

struct FastqRecord {
    read: Record,
    read2: Option<Record>,
    hash: u64,
}

impl FastqRecord {
    fn from_record(record: Record) -> Self {
        let hash = XxHash64::oneshot(0, record.sequence().as_bytes());
        Self {
            read: record,
            read2: None,
            hash,
        }
    }

    fn from_pair(read1: Record, read2: Record) -> Self {
        let s1 = read1.sequence();
        let s2 = read2.sequence();
        let hash = XxHash64::oneshot(0, &[s1, s2].concat());
        Self {
            read: read1,
            read2: Some(read2),
            hash,
        }
    }
}

// Maintains backward compatibility with the original API
pub struct DuplicateRemover<R>
where
    R: std::io::BufRead,
{
    fastq1: fastq::Reader<R>,
    fastq2: Option<fastq::Reader<R>>,
    seen: HashSet<u64>,
}

impl<R> DuplicateRemover<R>
where
    R: std::io::BufRead,
{
    fn new(fastq1: fastq::Reader<R>, fastq2: Option<fastq::Reader<R>>) -> Self {
        Self {
            fastq1,
            fastq2,
            seen: HashSet::new(),
        }
    }
}

impl DuplicateRemover<Box<dyn std::io::BufRead>> {
    pub fn from_fastq_paths<P>(fastq1: P, fastq2: Option<P>) -> Result<Self>
    where
        P: AsRef<Path> + Clone,
    {
        let fastq1 = crate::utils::get_fastq_reader(fastq1)?;
        let fastq2 = match fastq2 {
            Some(fastq2) => Some(crate::utils::get_fastq_reader(fastq2)?),
            None => None,
        };
        Ok(Self::new(fastq1, fastq2))
    }

    // Synchronous wrapper around the async function for compatibility
    pub fn deduplicate<P>(&mut self, output1: P, output2: Option<P>) -> Result<FastqDeduplicationStats>
    where
        P: AsRef<Path> + Clone + Send + 'static,
    {
        let rt = tokio::runtime::Runtime::new()?;
        rt.block_on(async {
            // Convert sync inputs to async inputs
            let async_result = AsyncDuplicateRemover::from_fastq_paths(
                self.fastq1

            ).await?;
            
            // Use the shared HashSet
            let mut async_remover = async_result;
            async_remover.seen = self.seen.clone();
            
            let stats = if self.fastq2.is_some() {
                async_remover.deduplicate_paired().await?
            } else {
                async_remover.deduplicate_single().await?
            };
            
            // Update the original set with newly seen hashes
            self.seen = async_remover.seen;
            
            Ok(stats)
        })
    }
}

// New async implementation
pub struct AsyncDuplicateRemover {
    fastq1_path: std::path::PathBuf,
    fastq2_path: Option<std::path::PathBuf>,
    output1_path: std::path::PathBuf,
    output2_path: Option<std::path::PathBuf>,
    seen: HashSet<u64>,
}

impl AsyncDuplicateRemover {
    pub async fn from_fastq_paths<P>(
        fastq1: P, 
        fastq2: Option<P>,
        output1: P,
        output2: Option<P>
    ) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(Self {
            fastq1_path: fastq1.as_ref().to_path_buf(),
            fastq2_path: fastq2.map(|p| p.as_ref().to_path_buf()),
            output1_path: output1.as_ref().to_path_buf(),
            output2_path: output2.map(|p| p.as_ref().to_path_buf()),
            seen: HashSet::new(),
        })
    }

    pub async fn deduplicate(&mut self) -> Result<FastqDeduplicationStats> {
        if self.fastq2_path.is_some() {
            self.deduplicate_paired().await
        } else {
            self.deduplicate_single().await
        }
    }

    async fn deduplicate_single(&mut self) -> Result<FastqDeduplicationStats> {
        // Open the input file with tokio
        let file = File::open(&self.fastq1_path).await?;
        let reader = BufReader::new(file);
        
        // Create noodles async reader
        let mut fastq_reader = AsyncReader::new(reader);
        
        // Open output file with tokio
        let output_file = File::create(&self.output1_path).await?;
        let writer = Arc::new(TokioMutex::new(output_file));
        
        let seen = Arc::new(TokioMutex::new(&mut self.seen));
        let stats = Arc::new(TokioMutex::new(FastqDeduplicationStats::new()));
        
        // Process records in parallel batches
        let mut tasks = Vec::new();
        let mut batch = Vec::with_capacity(BATCH_SIZE);
        
        let mut records_stream = fastq_reader.records();
        
        while let Some(record_result) = records_stream.next().await {
            let record = record_result?;
            batch.push(record);
            
            if batch.len() >= BATCH_SIZE {
                let current_batch = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE));
                let writer_clone = Arc::clone(&writer);
                let seen_clone = Arc::clone(&seen);
                let stats_clone = Arc::clone(&stats);
                
                let task = tokio::spawn(async move {
                    process_batch_single(current_batch, writer_clone, seen_clone, stats_clone).await
                });
                
                tasks.push(task);
            }
        }
        
        // Process any remaining records
        if !batch.is_empty() {
            let writer_clone = Arc::clone(&writer);
            let seen_clone = Arc::clone(&seen);
            let stats_clone = Arc::clone(&stats);
            
            let task = tokio::spawn(async move {
                process_batch_single(batch, writer_clone, seen_clone, stats_clone).await
            });
            
            tasks.push(task);
        }
        
        // Wait for all tasks to complete
        for task in tasks {
            task.await??;
        }
        
        // Extract final stats
        let final_stats = Arc::try_unwrap(stats)
            .expect("Failed to unwrap stats Arc")
            .into_inner()
            .expect("Failed to unwrap Mutex");
        
        Ok(final_stats)
    }

    async fn deduplicate_paired(&mut self) -> Result<FastqDeduplicationStats> {
        // Open the input files with tokio
        let file1 = File::open(&self.fastq1_path).await?;
        let reader1 = BufReader::new(file1);
        let mut fastq_reader1 = AsyncReader::new(reader1);
        
        let file2 = File::open(self.fastq2_path.as_ref().unwrap()).await?;
        let reader2 = BufReader::new(file2);
        let mut fastq_reader2 = AsyncReader::new(reader2);
        
        // Open output files with tokio
        let output_file1 = File::create(&self.output1_path).await?;
        let writer1 = Arc::new(TokioMutex::new(output_file1));
        
        let output_file2 = File::create(self.output2_path.as_ref().unwrap()).await?;
        let writer2 = Arc::new(TokioMutex::new(output_file2));
        
        let seen = Arc::new(TokioMutex::new(&mut self.seen));
        let stats = Arc::new(TokioMutex::new(FastqDeduplicationStats::new()));
        
        // Process records in parallel batches
        let mut tasks = Vec::new();
        let mut batch1 = Vec::with_capacity(BATCH_SIZE);
        let mut batch2 = Vec::with_capacity(BATCH_SIZE);
        
        let mut records_stream1 = fastq_reader1.records();
        let mut records_stream2 = fastq_reader2.records();
        
        // Process paired records in batches
        let mut i = 0;
        loop {
            let record1 = match records_stream1.next().await {
                Some(r) => r?,
                None => break,
            };
            
            let record2 = match records_stream2.next().await {
                Some(r) => r?,
                None => {
                    return Err(anyhow::anyhow!("Mismatched number of records in paired files"));
                }
            };
            
            batch1.push(record1);
            batch2.push(record2);
            i += 1;
            
            if i >= BATCH_SIZE {
                let current_batch1 = std::mem::replace(&mut batch1, Vec::with_capacity(BATCH_SIZE));
                let current_batch2 = std::mem::replace(&mut batch2, Vec::with_capacity(BATCH_SIZE));
                
                let writer1_clone = Arc::clone(&writer1);
                let writer2_clone = Arc::clone(&writer2);
                let seen_clone = Arc::clone(&seen);
                let stats_clone = Arc::clone(&stats);
                
                let task = tokio::spawn(async move {
                    process_batch_paired(
                        current_batch1, 
                        current_batch2, 
                        writer1_clone, 
                        writer2_clone, 
                        seen_clone, 
                        stats_clone
                    ).await
                });
                
                tasks.push(task);
                i = 0;
            }
        }
        
        // Process any remaining records
        if !batch1.is_empty() {
            let writer1_clone = Arc::clone(&writer1);
            let writer2_clone = Arc::clone(&writer2);
            let seen_clone = Arc::clone(&seen);
            let stats_clone = Arc::clone(&stats);
            
            let task = tokio::spawn(async move {
                process_batch_paired(
                    batch1, 
                    batch2, 
                    writer1_clone, 
                    writer2_clone, 
                    seen_clone, 
                    stats_clone
                ).await
            });
            
            tasks.push(task);
        }
        
        // Wait for all tasks to complete
        for task in tasks {
            task.await??;
        }
        
        // Extract final stats
        let final_stats = Arc::try_unwrap(stats)
            .expect("Failed to unwrap stats Arc")
            .into_inner()
            ;
        
        Ok(final_stats)
    }
}

// Helper function to process a batch of single-end records
async fn process_batch_single(
    batch: Vec<Record>,
    writer: Arc<TokioMutex<File>>,
    seen: Arc<TokioMutex<&mut HashSet<u64>>>,
    stats: Arc<TokioMutex<FastqDeduplicationStats>>
) -> Result<()> {
    let mut local_stats = FastqDeduplicationStats::new();
    let mut write_buffer = Vec::new();
    let mut new_hashes = Vec::new();
    
    // Process each record locally first to minimize mutex lock time
    for record in batch {
        local_stats.increment_total();
        
        let fastq_record = FastqRecord::from_record(record);
        let hash = fastq_record.hash;
        
        // Lock the shared hashset to check if we've seen this hash
        let is_seen = {
            let seen_lock = seen.lock().await;
            seen_lock.contains(&hash)
        };
        
        if !is_seen {
            local_stats.increment_unique();
            new_hashes.push(hash);
            
            // Format the record for writing
            let mut buf = Vec::new();
            write_fastq_record(&fastq_record.read, &mut buf)?;
            write_buffer.extend(buf);
        } else {
            local_stats.increment_duplicate();
        }
    }
    
    // Update the shared hashset with new hashes
    {
        let mut seen_lock = seen.lock().await;
        for hash in new_hashes {
            seen_lock.insert(hash);
        }
    }
    
    // Write all unique records
    if !write_buffer.is_empty() {
        let mut writer_lock = writer.lock().await;
        writer_lock.write_all(&write_buffer).await?;
    }
    
    // Update stats
    {
        let mut stats_lock = stats.lock().await;
        stats_lock.merge(&local_stats);
    }
    
    Ok(())
}

// Helper function to process a batch of paired-end records
async fn process_batch_paired(
    batch1: Vec<Record>,
    batch2: Vec<Record>,
    writer1: Arc<TokioMutex<File>>,
    writer2: Arc<TokioMutex<File>>,
    seen: Arc<TokioMutex<&mut HashSet<u64>>>,
    stats: Arc<TokioMutex<FastqDeduplicationStats>>
) -> Result<()> {
    let mut local_stats = FastqDeduplicationStats::new();
    let mut write_buffer1 = Vec::new();
    let mut write_buffer2 = Vec::new();
    let mut new_hashes = Vec::new();
    
    // Process each pair locally first to minimize mutex lock time
    for (record1, record2) in batch1.into_iter().zip(batch2.into_iter()) {
        local_stats.increment_total();
        
        let fastq_record = FastqRecord::from_pair(record1, record2);
        let hash = fastq_record.hash;
        
        // Lock the shared hashset to check if we've seen this hash
        let is_seen = {
            let seen_lock = seen.lock().await;
            seen_lock.contains(&hash)
        };
        
        if !is_seen {
            local_stats.increment_unique();
            new_hashes.push(hash);
            
            // Format the records for writing
            let mut buf1 = Vec::new();
            write_fastq_record(&fastq_record.read, &mut buf1)?;
            write_buffer1.extend(buf1);
            
            if let Some(read2) = fastq_record.read2 {
                let mut buf2 = Vec::new();
                write_fastq_record(&read2, &mut buf2)?;
                write_buffer2.extend(buf2);
            }
        } else {
            local_stats.increment_duplicate();
        }
    }
    
    // Update the shared hashset with new hashes
    {
        let mut seen_lock = seen.lock().await;
        for hash in new_hashes {
            seen_lock.insert(hash);
        }
    }
    
    // Write all unique records (first file)
    if !write_buffer1.is_empty() {
        let mut writer_lock = writer1.lock().await;
        writer_lock.write_all(&write_buffer1).await?;
    }
    
    // Write all unique records (second file)
    if !write_buffer2.is_empty() {
        let mut writer_lock = writer2.lock().await;
        writer_lock.write_all(&write_buffer2).await?;
    }
    
    // Update stats
    {
        let mut stats_lock = stats.lock().await;
        stats_lock.merge(&local_stats);
    }
    
    Ok(())
}

// Helper function to format a FASTQ record to bytes
fn write_fastq_record(record: &Record, buf: &mut Vec<u8>) -> Result<()> {
    writeln!(buf, "@{}", String::from_utf8_lossy(record.name()))?;
    writeln!(buf, "{}", String::from_utf8_lossy(record.sequence()))?;
    writeln!(buf, "+")?;
    writeln!(buf, "{}", String::from_utf8_lossy(record.quality_scores()))?;
    Ok(())
}

// Function to process multiple FASTQ files in parallel
pub async fn process_multiple_files<P>(
    input_files: Vec<(P, Option<P>)>,
    output_files: Vec<(P, Option<P>)>
) -> Result<Vec<FastqDeduplicationStats>>
where
    P: AsRef<Path> + Clone + Send + 'static,
{
    let mut tasks = Vec::new();
    
    for ((input1, input2), (output1, output2)) in input_files.into_iter().zip(output_files.into_iter()) {
        let task = task::spawn(async move {
            let mut deduplicator = AsyncDuplicateRemover::from_fastq_paths(
                input1, input2, output1, output2
            ).await?;
            deduplicator.deduplicate().await
        });
        tasks.push(task);
    }
    
    let mut results = Vec::new();
    for task in tasks {
        results.push(task.await??);
    }
    
    Ok(results)
}