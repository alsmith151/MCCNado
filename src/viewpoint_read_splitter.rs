use anyhow::Result;
use bio::io::fasta::Sequence;
use bio::utils;
use bstr::ByteSlice;
use flate2;
use log::{info, warn};
use noodles::fastq;
use noodles::fastq::record::Definition;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::header::record::value::map::header;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::de;
use std::collections::HashSet;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use bamnado::BamStats;

use crate::utils::{
    get_fastq_reader, get_fastq_writer, FlashedStatus, ReadNumber, Segment,
    SegmentMetadata, Strand, ViewpointPosition, SegmentType, SegmentPositions,
};

struct ViewpointRead<'a> {
    viewpoint: &'a str,
    read: noodles::bam::Record,
    flashed_status: FlashedStatus,
    minimum_segment_length: usize,
}

impl<'a> ViewpointRead<'a> {
    fn new(
        viewpoint_name: &'a str,
        read: noodles::bam::Record,
        flashed_status: FlashedStatus,
        minimum_segment_length: usize,
    ) -> Self {
        Self {
            viewpoint: viewpoint_name,
            read,
            flashed_status,
            minimum_segment_length,
        }
    }

    fn name(&self) -> Result<String> {
        match self.read.name() {
            Some(name) => Ok(name.to_string()),
            None => Err(anyhow::anyhow!("Read name not found")),
        }
    }

    fn is_viewpoint_read(&self) -> bool {
        !self.read.flags().is_unmapped()
    }

    fn strand(&self) -> Strand {
        if self.read.flags().is_reverse_complemented() {
            Strand::NEGATIVE
        } else {
            Strand::POSITIVE
        }
    }

    fn viewpoint(&self) -> Option<String> {
        if !self.is_viewpoint_read() {
            return None;
        }

        match self.read.reference_sequence_id() {
            Some(reference_sequence_id) => match reference_sequence_id {
                Ok(rid) => Some(rid.to_string()),
                Err(_) => {
                    return None;
                }
            },
            None => None,
        }
    }

    fn read_number(&self) -> ReadNumber {
        match self.flashed_status {
            FlashedStatus::FLASHED => ReadNumber::FLASHED,
            FlashedStatus::UNFLASHED => {
                if self.read.flags().is_first_segment() {
                    ReadNumber::ONE
                } else {
                    ReadNumber::TWO
                }
            }
        }
    }

    fn segments(&self) -> Option<Result<Vec<Segment<fastq::Record>>>> {
        if !self.is_viewpoint_read() {
            return None;
        }

        let mut segment_positions = SegmentPositions::default();

        // Get the CIGAR string
        let cigar = self.read.cigar();
        let mut start = 0;
        // let mut is_matched = false;

        let mut current_segment = SegmentType::LEFT;

        for op in cigar.iter() {
            let op = match op {
                Ok(op) => op,
                Err(_) => return None,
            };

            // Until a match is found, we are in the left segment
            // If bases have been soft-clipped, they are part of the left segment
            // Deletions == no change in position, insertions increase position to end of insertion
            // When a match is found, we are in the viewpoint segment
            // Same rules apply for Insertions and Deletions
            // When soft-clipped bases are found, we are in the right segment

            if op.kind() == Kind::SoftClip && current_segment == SegmentType::LEFT {
                segment_positions.set_left((start, start + op.len() as usize));
                start += op.len() as usize;
            } else if op.kind() == Kind::Match && current_segment == SegmentType::LEFT {
                segment_positions.set_viewpoint((start, start + op.len() as usize));
                current_segment = SegmentType::VIEWPOINT;
                start += op.len() as usize;
            } else if op.kind() == Kind::SoftClip && current_segment == SegmentType::VIEWPOINT {
                segment_positions.set_right((start, start + op.len() as usize));
                current_segment = SegmentType::RIGHT;
                start += op.len() as usize;
            } else if op.kind() == Kind::Insertion {
                start += op.len() as usize;
            } else if op.kind() == Kind::Match {
                let seg = match current_segment {
                    SegmentType::LEFT => &mut segment_positions.left(),
                    SegmentType::VIEWPOINT => &mut segment_positions.viewpoint(),
                    SegmentType::RIGHT => &mut segment_positions.right(),
                };
                seg.1 += op.len() as usize;
            }
        }

        let mut segments = Vec::new();
        let sequence = self.read.sequence();
        let quality_scores = self.read.quality_scores();

        for (segment_type, positions) in segment_positions.into_iter() {
            let (start, end) = positions;

            let end = match end > sequence.len() {
                true => {
                    warn!(
                        "Segment end is greater than sequence length for read {}",
                        self.name().unwrap()
                    );
                    sequence.len()
                }
                false => end,
            };

            let sequence: Vec<u8> = sequence
                .iter()
                .enumerate()
                .skip(start)
                .take(end - start)
                .map(|(i, b)| sequence.get(i).unwrap())
                .collect();

            let quality_scores: &[u8] = &quality_scores.as_ref()[start..end]
                .iter()
                .map(|&n| n + 33)
                .collect::<Vec<u8>>();

            let metadata = SegmentMetadata::from_parts(
                self.name().unwrap().as_str(),
                &self.viewpoint,
                ViewpointPosition::from_segment_type(segment_type),
                self.read_number(),
                self.flashed_status,
            );

            if sequence.len() < self.minimum_segment_length {
                continue;
            } else {
                segments.push(Segment::<fastq::Record>::from_metadata(metadata, &sequence, quality_scores));
            }
        }

        Some(Ok(segments))
    }
}

pub struct ReadSplitterOptions {
    flashed_status: FlashedStatus,
    minimum_segment_length: usize,
}

impl Default for ReadSplitterOptions {
    fn default() -> Self {
        Self {
            flashed_status: FlashedStatus::FLASHED,
            minimum_segment_length: 18,
        }
    }
}

pub struct ReadSplitter {
    bam_path: PathBuf,
    options: ReadSplitterOptions,
}

impl ReadSplitter {
    pub fn new(bam_path: &str, options: ReadSplitterOptions) -> Self {
        Self {
            bam_path: PathBuf::from(bam_path),
            options,
        }
    }

    fn header(&self) -> Result<noodles::sam::Header> {
        let header_samtools = std::process::Command::new("samtools")
            .arg("view")
            .arg("-H")
            .arg(self.bam_path.clone())
            .output()
            .expect("Failed to run samtools")
            .stdout;

        let header_str =
            String::from_utf8(header_samtools).expect("Failed to convert header to string");

        // Slight hack here for CellRanger BAM files that are missing the version info
        let header_string =
            header_str.replace("@HD\tSO:coordinate\n", "@HD\tVN:1.6\tSO:coordinate\n");
        let header_str = header_string.as_bytes();
        let mut reader = noodles::sam::io::Reader::new(header_str);
        let header = reader
            .read_header()
            .expect("Failed to read header with samtools");
        Ok(header)
    }

    pub fn split_reads(&self, outfile: &str) -> Result<()> {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(self.bam_path.clone()).expect("Failed to build indexed reader");
        let mut writer = get_fastq_writer(outfile).expect("Failed to create fastq writer");

        let header = self.header()?;
        // Get the chromosome sizes
        let chromsizes = header
            .reference_sequences()
            .iter()
            .map(|(name, seq)| (name.to_string(), seq.length().get() as u64))
            .collect::<std::collections::HashMap<_, _>>();

        let query_regions = chromsizes.iter().map(|(name, size)| {
            let start = noodles::core::Position::try_from(1).unwrap();
            let end = noodles::core::Position::try_from(*size as usize).unwrap();
            noodles::core::Region::new(name.to_string(), start..=end)
        });

        let mut counter = 0;
        for region in query_regions {
            let mut records = reader.query(&header, &region).expect("Failed to query region");
            while let Some(record) = records.next() {
                counter += 1;
                if counter % 100_000 == 0 {
                    info!("Processed {} reads", counter);
                }

                let viewpoint_read = ViewpointRead::new(
                    region.name().to_str().unwrap(),
                    record?,
                    self.options.flashed_status,
                    self.options.minimum_segment_length,
                );
                if let Some(segments) = viewpoint_read.segments() {
                    let segments = segments?;
                    for segment in segments {
                        writer.write_record(segment.record())?;
                    }
                }
            }
        }


        Ok(())
    }
}
