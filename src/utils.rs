use anyhow::{Context, Result};
use bstr::ByteSlice;
use noodles::fastq;
use noodles::fastq::record::Definition;
use std::path::Path;

use noodles::bam;

pub fn get_fastq_reader<P>(
    fname: P,
) -> Result<noodles::fastq::io::Reader<Box<dyn std::io::BufRead>>>
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
    Flashed = 1,
    #[allow(dead_code)]
    Unflashed = 0,
}

impl std::fmt::Display for FlashedStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            FlashedStatus::Flashed => write!(f, "1"),
            FlashedStatus::Unflashed => write!(f, "0"),
        }
    }
}

impl FlashedStatus {
    #[allow(dead_code)]
    pub fn from_str(s: &str) -> Self {
        match s {
            "1" => FlashedStatus::Flashed,
            "0" => FlashedStatus::Unflashed,
            _ => panic!("Invalid flashed status"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReadNumber {
    One = 1,
    Two = 2,
    Flashed = 3,
}

impl std::fmt::Display for ReadNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ReadNumber::One => write!(f, "1"),
            ReadNumber::Two => write!(f, "2"),
            ReadNumber::Flashed => write!(f, "3"),
        }
    }
}

impl ReadNumber {
    #[allow(dead_code)]
    fn from_str(s: &str) -> Self {
        match s {
            "1" => ReadNumber::One,
            "2" => ReadNumber::Two,
            _ => panic!("Invalid read number"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Strand {
    Positive = 1,
    Negative = -1,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::Positive => write!(f, "1"),
            Strand::Negative => write!(f, "-1"),
        }
    }
}

impl Strand {
    #[allow(dead_code)]
    fn from_str(s: &str) -> Self {
        match s {
            "1" => Strand::Positive,
            "-1" => Strand::Negative,
            _ => panic!("Invalid strand"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SegmentType {
    Left,
    Viewpoint,
    Right,
}

impl std::fmt::Display for SegmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SegmentType::Left => write!(f, "left"),
            SegmentType::Viewpoint => write!(f, "viewpoint"),
            SegmentType::Right => write!(f, "right"),
        }
    }
}

impl SegmentType {
    #[allow(dead_code)]
    fn from_str(s: &str) -> Self {
        match s {
            "left" => SegmentType::Left,
            "viewpoint" => SegmentType::Viewpoint,
            "right" => SegmentType::Right,
            _ => panic!("Invalid segment type"),
        }
    }

    pub fn from_viewpoint_position(viewpoint_position: ViewpointPosition) -> Self {
        match viewpoint_position {
            ViewpointPosition::Start => SegmentType::Right,
            ViewpointPosition::End => SegmentType::Left,
            ViewpointPosition::All => SegmentType::Viewpoint,
            ViewpointPosition::None => panic!("Invalid viewpoint position"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ViewpointPosition {
    Start = 5,
    End = 3,
    All = 1,
    None = 0,
}

impl std::fmt::Display for ViewpointPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ViewpointPosition::Start => write!(f, "start"),
            ViewpointPosition::End => write!(f, "end"),
            ViewpointPosition::All => write!(f, "all"),
            ViewpointPosition::None => write!(f, "none"),
        }
    }
}

impl ViewpointPosition {
    fn from_str(s: &str) -> Self {
        match s {
            "start" => ViewpointPosition::Start,
            "end" => ViewpointPosition::End,
            "all" => ViewpointPosition::All,
            "none" => ViewpointPosition::None,
            _ => panic!("Invalid viewpoint position"),
        }
    }
}

impl ViewpointPosition {
    pub fn from_segment_type(segment_type: SegmentType) -> Self {
        match segment_type {
            SegmentType::Left => ViewpointPosition::End,
            SegmentType::Viewpoint => ViewpointPosition::All,
            SegmentType::Right => ViewpointPosition::Start,
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
        self.viewpoint()
            .split_once("-")
            .context("No viewpoint coordinate")
            .expect("Error splitting oligo coords")
            .1
    }

    pub fn viewpoint_name(&self) -> &str {
        self.viewpoint()
            .split_once("-")
            .context("No viewpoint coordinate")
            .expect("Error splitting oligo coords")
            .0
    }

    pub fn viewpoint_position(&self) -> ViewpointPosition {
        ViewpointPosition::from_str(self.name.split("__").nth(2).unwrap())
    }

    #[allow(dead_code)]
    pub fn read_number(&self) -> ReadNumber {
        ReadNumber::from_str(self.name.split("__").nth(3).unwrap())
    }

    #[allow(dead_code)]
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
pub struct Segment<R> {
    #[allow(dead_code)]
    metadata: SegmentMetadata,
    record: R,
}

impl<R> Segment<R> {
    #[allow(dead_code)]
    fn new(metadata: SegmentMetadata, record: R) -> Self {
        Self { metadata, record }
    }

    #[allow(dead_code)]
    pub fn metadata(&self) -> &SegmentMetadata {
        &self.metadata
    }

    pub fn record(&self) -> &R {
        &self.record
    }
}

impl Segment<fastq::Record> {
    pub fn from_metadata(
        metadata: SegmentMetadata,
        sequence: &[u8],
        quality_scores: &[u8],
    ) -> Self {
        let name = metadata.name.as_bytes();
        let record = fastq::Record::new(Definition::new(name, ""), sequence, quality_scores);
        Self { metadata, record }
    }
}

impl Segment<bam::Record> {
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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

    #[allow(dead_code)]
    pub fn set_positions(
        &mut self,
        viewpoint: (usize, usize),
        left: (usize, usize),
        right: (usize, usize),
    ) {
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
            Some((SegmentType::Left, self.left))
        } else if self.current_pos == 1 {
            self.current_pos += 1;
            Some((SegmentType::Viewpoint, self.viewpoint))
        } else if self.current_pos == 2 {
            self.current_pos += 1;
            Some((SegmentType::Right, self.right))
        } else {
            None
        }
    }
}
