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

pub fn get_fastq_reader<P>(fname: P) -> Result<noodles::fastq::io::Reader<Box<dyn std::io::BufRead>>>
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

    Ok(noodles::fastq::io::Reader::new(buffer))
}

pub fn get_fastq_writer<P>(fname: P) -> Result<noodles::fastq::io::Writer<Box<dyn std::io::Write>>>
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

    pub fn viewpoint_name(&self) -> &str {
        self.viewpoint().split_once("-").context("No viewpoint coordinate").expect("Error splitting oligo coords").0
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
