use std::path::Path;
use std::io::{Write};
use std::collections::HashSet;
use anyhow::Result;
use noodles::fastq;
use twox_hash::XxHash64;
use bstr::ByteSlice;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::utils::{get_fastq_reader, get_fastq_writer};


#[derive(Debug, serde::Serialize, serde::Deserialize, Clone, IntoPyObject)]
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

    fn print(&self) {
        println!("Total reads: {}", self.total_reads);
        println!("Unique reads: {}", self.unique_reads);
        println!("Duplicate reads: {}", self.duplicate_reads);
    }
}

// impl IntoPy<PyObject> for FastqDeduplicationStats {
//     fn into_py(self, py: Python) -> PyObject {
//         let dict = PyDict::new_bound(py);
//         dict.set_item("total_reads", self.total_reads).unwrap();
//         dict.set_item("unique_reads", self.unique_reads).unwrap();
//         dict.set_item("duplicate_reads", self.duplicate_reads).unwrap();
//         dict.into()
//     }
// }


struct FastqRecord {
    read: fastq::Record,
    read2: Option<fastq::Record>,
    hash: u64,
}

impl FastqRecord {
    fn from_record(record: fastq::Record) -> Self {
        let hash = XxHash64::oneshot(0, record.sequence().as_bytes());
        Self {
            read: record,
            read2: None,
            hash,
        }
    }

    fn from_pair(read1: fastq::Record, read2: fastq::Record) -> Self {
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

pub struct DuplicateRemover<R>
where
    R: std::io::BufRead,
{
    fastq1: noodles::fastq::io::Reader<R>,
    fastq2: Option<noodles::fastq::io::Reader<R>>,
    seen: HashSet<u64>,
}

impl<R> DuplicateRemover<R>
where
    R: std::io::BufRead,
{
    fn new(fastq1: noodles::fastq::io::Reader<R>, fastq2: Option<noodles::fastq::io::Reader<R>>) -> Self {
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
        let fastq1 = get_fastq_reader(fastq1)?;
        let fastq2 = match fastq2 {
            Some(fastq2) => Some(get_fastq_reader(fastq2)?),
            None => None,
        };
        Ok(Self::new(fastq1, fastq2))
    }

    pub fn deduplicate<P>(&mut self, output1: P, output2: Option<P>) -> Result<FastqDeduplicationStats>
    where
        P: AsRef<Path> + Clone,
    {
        let mut writer1 = get_fastq_writer(output1)?;
        let mut writer2 = match output2 {
            Some(output2) => Some(get_fastq_writer(output2)?),
            None => None,
        };

        if let Some(ref mut fastq2_reader) = self.fastq2 {
            let stats = self.deduplicate_paired(&mut writer1, writer2.as_mut().expect("Paired-end requires two output files"))?;
            return Ok(stats)
        } else {
            let stats = self.deduplicate_single(&mut writer1)?;
            return Ok(stats)
        }
    }

    fn deduplicate_single(&mut self, writer: &mut fastq::io::Writer<Box<dyn Write>>) -> Result<FastqDeduplicationStats> {
        
        let mut stats = FastqDeduplicationStats::new();
        for (ii, record) in self.fastq1.records().enumerate() {

            let record = record?;
            stats.increment_total();

            let fastq_record = FastqRecord::from_record(record);
            if self.seen.insert(fastq_record.hash) {
                stats.increment_unique();
                writer.write_record(&fastq_record.read)?;
            } else {
                stats.increment_duplicate();
            }
        }
        Ok(stats)
    }

    fn deduplicate_paired(
        &mut self,
        writer1: &mut fastq::io::Writer<Box<dyn Write>>,
        writer2: &mut fastq::io::Writer<Box<dyn Write>>,
    ) -> Result<FastqDeduplicationStats> {

        let mut stats = FastqDeduplicationStats::new();
        if let Some(ref mut fastq2_reader) = self.fastq2 {
            for (record1, record2) in self.fastq1.records().zip(fastq2_reader.records()) {
                let record1 = record1?;
                let record2 = record2?;
                stats.increment_total();
                let fastq_record = FastqRecord::from_pair(record1, record2);
                if self.seen.insert(fastq_record.hash) {
                    stats.increment_unique();
                    
                    writer1.write_record(&fastq_record.read)?;
                    if let Some(read2) = fastq_record.read2 {
                        writer2.write_record(&read2)?;
                    }
                } else {
                    stats.increment_duplicate();
                }
            }
        }
        Ok(stats)
    }
}
