use anyhow::{Context, Result};
use bio::bio_types::annot::contig::Contig;
use bio::bio_types::strand::ReqStrand;
use bstr::ByteSlice;
use flate2;
use log::{debug, info, warn};
use noodles::fastq;
use noodles::fastq::record::Definition;
use std::io::{BufRead, Write};
use std::path::Path;
use std::path::PathBuf;
use std::collections::HashMap;

use noodles::core::{Position, Region};
use noodles::{bam, sam};

pub fn get_fastq_reader<P>(fname: P) -> Result<fastq::Reader<Box<dyn std::io::BufRead>>>
where
    P: AsRef<Path> + Clone,
{
    let f = std::fs::File::open(fname.clone())?;

    let buffer: Box<dyn std::io::BufRead> = match fname.as_ref().extension() {
        Some(ext) if ext == "gz" => {
            let gz = flate2::read::MultiGzDecoder::new(f);
            Box::new(std::io::BufReader::new(gz))
        }
        _ => Box::new(std::io::BufReader::new(f)),
    };

    Ok(fastq::Reader::new(buffer))
}

pub fn get_fastq_writer<P>(fname: P) -> Result<fastq::io::Writer<Box<dyn std::io::Write>>>
where
    P: AsRef<Path> + Clone,
{
    let f = std::fs::File::create(fname.clone())?;

    let buffer_size = 16 * 1024 * 1024; // 16 MB
    let f = std::io::BufWriter::with_capacity(buffer_size, f);

    let buffer: Box<dyn std::io::Write> = match fname.as_ref().extension() {
        Some(ext) => {
            if ext == "gz" {
                let gz = flate2::write::GzEncoder::new(f, flate2::Compression::default());
                Box::new(gz)
            } else {
                Box::new(f)
            }
        }
        None => Box::new(f),
    };
    Ok(fastq::io::Writer::new(buffer))
}



#[derive(Debug)]
struct ChromosomeStats {
    chrom: String,
    length: u64,
    mapped: u64,
    unmapped: u64,
}

pub struct BamStats {
    // BAM file stats

    // File path
    file_path: PathBuf,

    // Header
    header: sam::Header,

    // Contigs present in the BAM file
    contigs: Vec<Contig<String, ReqStrand>>,

    // Chromosome stats
    chrom_stats: HashMap<String, ChromosomeStats>,

    // Total number of reads
    n_reads: u64,
    // Number of reads mapped to the genome
    n_mapped: u64,
    // Number of reads unmapped to the genome
    n_unmapped: u64,
}


pub fn bam_header(file_path: PathBuf) -> Result<sam::Header> {
    // Read the header of a BAM file
    // If the noodles crate fails to read the header, we fall back to samtools
    // This is a bit of a hack but it works for now
    // The noodles crate is more strict about the header format than samtools

    // Check that the file exists
    if !file_path.exists() {
        return Err(anyhow::Error::from(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("File not found: {}", file_path.display()),
        )));
    };

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(file_path.clone())
        .expect("Failed to open file");

    let header = match reader.read_header() {
        std::result::Result::Ok(header) => header,
        Err(e) => {
            debug!(
                "Failed to read header using noodels falling back to samtools: {}",
                e
            );

            let header_samtools = std::process::Command::new("samtools")
                .arg("view")
                .arg("-H")
                .arg(file_path.clone())
                .output()
                .expect("Failed to run samtools")
                .stdout;

            let header_str =
                String::from_utf8(header_samtools).expect("Failed to convert header to string");

            // Slight hack here for CellRanger BAM files that are missing the version info
            let header_string =
                header_str.replace("@HD\tSO:coordinate\n", "@HD\tVN:1.6\tSO:coordinate\n");
            let header_str = header_string.as_bytes();
            let mut reader = sam::io::Reader::new(header_str);
            let header = reader
                .read_header()
                .expect("Failed to read header with samtools");
            header
        }
    };
    Ok(header)
}


impl BamStats {
    pub fn new(file_path: PathBuf) -> Result<Self> {
        // Read the header of the BAM file
        let header = bam_header(file_path.clone())?;

        // Get the contigs from the header
        let contigs = header
            .reference_sequences()
            .iter()
            .map(|(name, map)| {
                let name = name.to_string();
                let length = map.length().get();
                let contig = Contig::new(name, 1, length, ReqStrand::Forward);
                contig
            })
            .collect();

        // Get the index of the BAM file
        let bam_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(file_path.clone())
            .expect("Failed to open file");
        let index = bam_reader.index();

        // Get the chromosome stats
        let mut chrom_stats = HashMap::default();

        for ((reference_sequence_name_buf, reference_sequence), index_reference_sequence) in header
            .reference_sequences()
            .iter()
            .zip(index.reference_sequences())
        {
            let reference_sequence_name = std::str::from_utf8(reference_sequence_name_buf)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

            let (mapped_record_count, unmapped_record_count) = index_reference_sequence
                .metadata()
                .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
                .unwrap_or_default();

            let stats = ChromosomeStats {
                chrom: reference_sequence_name.to_string(),
                length: usize::from(reference_sequence.length()) as u64,
                mapped: mapped_record_count,
                unmapped: unmapped_record_count,
            };

            chrom_stats.insert(reference_sequence_name.to_string(), stats);
        }

        let mut unmapped_record_count = index.unplaced_unmapped_record_count().unwrap_or_default();

        let mut mapped_record_count = 0;

        for (_, stats) in chrom_stats.iter() {
            mapped_record_count += stats.mapped;
            unmapped_record_count += stats.unmapped;
        }

        let n_reads = mapped_record_count + unmapped_record_count;

        Ok(Self {
            file_path,
            header,
            contigs,
            chrom_stats,
            n_reads,
            n_mapped: mapped_record_count,
            n_unmapped: unmapped_record_count,
        })
    }

    pub fn estimate_genome_chunk_length(&self, bin_size: u64) -> Result<u64> {
        // genomeLength = sum(bamHandles[0].lengths)
        // max_reads_per_bp = max([float(x) / genomeLength for x in mappedList])
        // genomeChunkLength = int(min(5e6, int(2e6 / (max_reads_per_bp * len(bamHandles)))))
        // genomeChunkLength -= genomeChunkLength % tile_size

        let stats = &self.chrom_stats;
        let genome_length = stats.values().map(|x| x.length).sum::<u64>();
        let max_reads_per_bp = self.n_mapped as f64 / genome_length as f64;
        let genome_chunk_length = f64::from(5e6).min(2e6 / (max_reads_per_bp));

        let correction = genome_chunk_length % bin_size as f64;
        let genome_chunk_length = genome_chunk_length - correction;

        Ok(genome_chunk_length as u64)
    }

    pub fn genome_chunks(&self, bin_size: u64) -> Result<Vec<Region>> {
        let genome_chunk_length = self.estimate_genome_chunk_length(bin_size)?;
        let chrom_chunks = self
            .chrom_stats
            .iter()
            .map(|(chrom, stats)| {
                let mut chunks = Vec::new();
                let mut chunk_start = 1;
                let chrom_end = stats.length;

                while chunk_start <= chrom_end {
                    // Corrected to include the last position in the range
                    let chunk_end = chunk_start + genome_chunk_length - 1; // Adjust to ensure the chunk covers exactly genome_chunk_length positions
                    let chunk_end = chunk_end.min(chrom_end); // Ensure we do not exceed the chromosome length

                    let start = Position::try_from(chunk_start as usize).unwrap();
                    let end = Position::try_from(chunk_end as usize).unwrap();

                    let region = Region::new(&*chrom.clone(), start..=end);
                    chunks.push(region);
                    chunk_start = chunk_end + 1; // Corrected to start the next chunk right after the current chunk ends
                }
                chunks
            })
            .flatten()
            .collect();

        Ok(chrom_chunks)
    }

    pub fn ref_id_mapping(&self) -> HashMap<usize, String> {
        let mut ref_id_mapping = HashMap::default();
        for (i, (name, _)) in self.header.reference_sequences().iter().enumerate() {
            let name = std::str::from_utf8(name)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
                .unwrap()
                .to_string();
            ref_id_mapping.insert(i, name);
        }
        ref_id_mapping
    }

    pub fn chromsizes_ref_id(&self) -> Result<HashMap<usize, u64>> {
        let mut chromsizes = HashMap::default();
        for (i, (_, map)) in self.header.reference_sequences().iter().enumerate() {
            let length = map.length().get();
            chromsizes.insert(i, length as u64);
        }
        Ok(chromsizes)
    }

    pub fn chromsizes_ref_name(&self) -> Result<HashMap<String, u64>> {
        let mut chromsizes = HashMap::default();
        for (name, map) in self.header.reference_sequences() {
            let length = map.length().get();
            let name = std::str::from_utf8(name)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
                .unwrap()
                .to_string();
            chromsizes.insert(name, length as u64);
        }
        Ok(chromsizes)
    }

    pub fn write_chromsizes(&self, outfile: PathBuf) -> Result<()> {
        // let chrom_sizes: HashMap<String, u64> = self
        //     .chrom_stats
        //     .iter()
        //     .map(|(name, stats)| (name.to_owned(), stats.length.clone()))
        //     .collect();

        let chrom_sizes = self.chromsizes_ref_name()?;

        info!("Writing chromosome sizes to {}", outfile.display());
        let mut writer = std::io::BufWriter::new(std::fs::File::create(outfile)?);

        for (chrom, length) in chrom_sizes.iter() {
            writeln!(writer, "{}\t{}", chrom, length)?;
        }
        Ok(())
    }
}




#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FlashedStatus {
    FLASHED = 1,
    UNFLASHED = 0,
}

impl std::fmt::Display for FlashedStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            FlashedStatus::FLASHED => write!(f, "1"),
            FlashedStatus::UNFLASHED => write!(f, "0"),
        }
    }
}

impl FlashedStatus {
    pub fn from_str(s: &str) -> Self {
        match s {
            "1" => FlashedStatus::FLASHED,
            "0" => FlashedStatus::UNFLASHED,
            _ => panic!("Invalid flashed status"),
        }
    }
}



#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReadNumber {
    ONE = 1,
    TWO = 2,
    FLASHED = 3,
}

impl std::fmt::Display for ReadNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ReadNumber::ONE => write!(f, "1"),
            ReadNumber::TWO => write!(f, "2"),
            ReadNumber::FLASHED => write!(f, "3"),
        }
    }
}

impl ReadNumber {
    fn from_str(s: &str) -> Self {
        match s {
            "1" => ReadNumber::ONE,
            "2" => ReadNumber::TWO,
            _ => panic!("Invalid read number"),
        }
    }
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Strand {
    POSITIVE = 1,
    NEGATIVE = -1,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::POSITIVE => write!(f, "1"),
            Strand::NEGATIVE => write!(f, "-1"),
        }
    }
}

impl Strand {
    fn from_str(s: &str) -> Self {
        match s {
            "1" => Strand::POSITIVE,
            "-1" => Strand::NEGATIVE,
            _ => panic!("Invalid strand"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SegmentType {
    LEFT,
    VIEWPOINT,
    RIGHT,
}

impl std::fmt::Display for SegmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SegmentType::LEFT => write!(f, "left"),
            SegmentType::VIEWPOINT => write!(f, "viewpoint"),
            SegmentType::RIGHT => write!(f, "right"),
        }
    }
}

impl SegmentType {
    fn from_str(s: &str) -> Self {
        match s {
            "left" => SegmentType::LEFT,
            "viewpoint" => SegmentType::VIEWPOINT,
            "right" => SegmentType::RIGHT,
            _ => panic!("Invalid segment type"),
        }
    }

    pub fn from_viewpoint_position(viewpoint_position: ViewpointPosition) -> Self {
        match viewpoint_position {
            ViewpointPosition::START => SegmentType::RIGHT,
            ViewpointPosition::END => SegmentType::LEFT,
            ViewpointPosition::ALL => SegmentType::VIEWPOINT,
            ViewpointPosition::NONE => panic!("Invalid viewpoint position"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ViewpointPosition {
    START = 5,
    END = 3,
    ALL = 1,
    NONE = 0,
}

impl std::fmt::Display for ViewpointPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ViewpointPosition::START => write!(f, "start"),
            ViewpointPosition::END => write!(f, "end"),
            ViewpointPosition::ALL => write!(f, "all"),
            ViewpointPosition::NONE => write!(f, "none"),
        }
    }
}

impl ViewpointPosition {
    fn from_str(s: &str) -> Self {
        match s {
            "start" => ViewpointPosition::START,
            "end" => ViewpointPosition::END,
            "all" => ViewpointPosition::ALL,
            "none" => ViewpointPosition::NONE,
            _ => panic!("Invalid viewpoint position"),
        }
    }
}

impl ViewpointPosition {
    pub fn from_segment_type(segment_type: SegmentType) -> Self {
        match segment_type {
            SegmentType::LEFT => ViewpointPosition::END,
            SegmentType::VIEWPOINT => ViewpointPosition::ALL,
            SegmentType::RIGHT => ViewpointPosition::START,
        }
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct SegmentMetadata {
    name: String,
}

impl SegmentMetadata {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
        }
    }

    pub fn from_read_name(name: Option<&bstr::BStr>) -> Self {
        
        let name = match name {
            Some(name) => name,
            None => "UNKNOWN".into(),
        };

        Self {
            name: name.to_str().unwrap().to_string(),
        }
        
       
    }

    pub fn parent_id(&self) -> &str {
        self.name.split("__").next().unwrap()
    }

    pub fn viewpoint(&self) -> &str {
        self.name.split("__").nth(1).unwrap()
    }

    pub fn oligo_coordinates(&self) -> &str {
        self.viewpoint().split_once("-").context("No viewpoint coordinate").expect("Error splitting oligo coords").1
    }

    pub fn viewpoint_position(&self) -> ViewpointPosition {
        ViewpointPosition::from_str(self.name.split("__").nth(2).unwrap())
    }

    pub fn read_number(&self) -> ReadNumber {
        ReadNumber::from_str(self.name.split("__").nth(3).unwrap())
    }

    pub fn flashed_status(&self) -> FlashedStatus {
        FlashedStatus::from_str(self.name.split("__").nth(4).unwrap())
    }

    pub fn from_parts(
        parent_id: &str,
        viewpoint: &str,
        viewpoint_position: ViewpointPosition,
        read_number: ReadNumber,
        flashed_status: FlashedStatus,
    ) -> Self {
        Self {
            name: format!(
                "{}__{}__{}__{}__{}",
                parent_id, viewpoint, viewpoint_position, read_number, flashed_status
            ),
        }
    }
}

impl std::fmt::Display for SegmentMetadata {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

impl std::fmt::Debug for SegmentMetadata {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "ReporterName({})", self.name)
    }
}







#[derive(Clone, Debug)]
pub struct Segment<R>{
    metadata: SegmentMetadata,
    record: R,
}

impl <R>Segment<R> {
    fn new(metadata: SegmentMetadata, record: R) -> Self {
        Self { metadata, record }
    }

    pub fn metadata(&self) -> &SegmentMetadata {
        &self.metadata
    }

    pub fn record(&self) -> &R {
        &self.record
    }
}


impl Segment <fastq::Record> {
    pub fn from_metadata(metadata: SegmentMetadata, sequence: &[u8], quality_scores: &[u8]) -> Self {
        let name = metadata.name.as_bytes();
        let record = fastq::Record::new(Definition::new(name, ""), sequence, quality_scores);
        Self { metadata, record }
    }

}

impl Segment <bam::Record> {
    pub fn from_metadata(metadata: SegmentMetadata, record: bam::Record) -> Self {
        Self { metadata, record }
    }
}


#[derive(Debug)]
pub struct SegmentPositions {
    viewpoint: (usize, usize),
    left: (usize, usize),
    right: (usize, usize),

    current_pos: usize,
}

impl SegmentPositions {
    fn new(viewpoint: (usize, usize), left: (usize, usize), right: (usize, usize)) -> Self {
        Self {
            viewpoint,
            left,
            right,
            current_pos: 0,
        }
    }

    pub fn default() -> Self {
        Self {
            viewpoint: (0, 0),
            left: (0, 0),
            right: (0, 0),
            current_pos: 0,
        }
    }

    pub fn viewpoint(&self) -> (usize, usize) {
        self.viewpoint
    }

    pub fn left(&self) -> (usize, usize) {
        self.left
    }

    pub fn right(&self) -> (usize, usize) {
        self.right
    }

    pub fn set_viewpoint(&mut self, viewpoint: (usize, usize)) {
        self.viewpoint = viewpoint;
    }

    pub fn set_left(&mut self, left: (usize, usize)) {
        self.left = left;
    }

    pub fn set_right(&mut self, right: (usize, usize)) {
        self.right = right;
    }

    pub fn set_positions(&mut self, viewpoint: (usize, usize), left: (usize, usize), right: (usize, usize)) {
        self.viewpoint = viewpoint;
        self.left = left;
        self.right = right;
    }

}

impl Iterator for SegmentPositions {
    type Item = (SegmentType, (usize, usize));

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_pos == 0 {
            self.current_pos += 1;
            return Some((SegmentType::LEFT, self.left));
        } else if self.current_pos == 1 {
            self.current_pos += 1;
            return Some((SegmentType::VIEWPOINT, self.viewpoint));
        } else if self.current_pos == 2 {
            self.current_pos += 1;
            return Some((SegmentType::RIGHT, self.right));
        } else {
            return None;
        }
    }
}
