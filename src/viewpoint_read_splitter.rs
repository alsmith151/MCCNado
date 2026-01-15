use anyhow::Result;
use bstr::ByteSlice;
use log::info;
use noodles::fastq;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::PathBuf;

use crate::utils::{
    get_fastq_writer, FlashedStatus, ReadNumber, Segment, SegmentMetadata, SegmentPositions,
    SegmentType, Strand, ViewpointPosition,
};

pub struct ViewpointRead<'a> {
    pub viewpoint: &'a str,
    pub read: noodles::sam::alignment::RecordBuf,
    pub flashed_status: FlashedStatus,
    pub minimum_segment_length: usize,
}

impl<'a> ViewpointRead<'a> {
    fn new(
        viewpoint_name: &'a str,
        read: noodles::sam::alignment::RecordBuf,
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

    #[allow(dead_code)]
    fn strand(&self) -> Strand {
        if self.read.flags().is_reverse_complemented() {
            Strand::Negative
        } else {
            Strand::Positive
        }
    }

    #[allow(dead_code)]
    fn viewpoint(&self) -> Option<String> {
        if !self.is_viewpoint_read() {
            return None;
        }

        self.read.reference_sequence_id().map(|rid| rid.to_string())
    }

    fn read_number(&self) -> ReadNumber {
        match self.flashed_status {
            FlashedStatus::Flashed => ReadNumber::Flashed,
            FlashedStatus::Unflashed => {
                if self.read.flags().is_first_segment() {
                    ReadNumber::One
                } else {
                    ReadNumber::Two
                }
            }
        }
    }

    fn parse_segment_positions(&self) -> Option<SegmentPositions> {
        let mut positions = SegmentPositions::default();
        let mut current_segment = SegmentType::Left;
        let mut offset = 0;

        for op in self.read.cigar().as_ref().iter() {
            let len = op.len();

            match (op.kind(), current_segment) {
                (Kind::SoftClip, SegmentType::Left) => {
                    positions.set_left((offset, offset + len));
                    offset += len;
                }
                (Kind::Match, SegmentType::Left) => {
                    positions.set_viewpoint((offset, offset + len));
                    current_segment = SegmentType::Viewpoint;
                    offset += len;
                }
                (Kind::SoftClip, SegmentType::Viewpoint) => {
                    positions.set_right((offset, offset + len));
                    current_segment = SegmentType::Right;
                    offset += len;
                }
                (Kind::Insertion, _) => offset += len,
                (Kind::Match, _) => {
                    let seg_range = match current_segment {
                        SegmentType::Left => &mut positions.left(),
                        SegmentType::Viewpoint => &mut positions.viewpoint(),
                        SegmentType::Right => &mut positions.right(),
                    };
                    seg_range.1 += len;
                }
                _ => {}
            }
        }
        Some(positions)
    }

    fn extract_segment(
        &self,
        segment_type: SegmentType,
        range: (usize, usize),
    ) -> Option<Segment<fastq::Record>> {
        let (start, end) = range;
        let seq_len = self.read.sequence().len();

        if start >= end || start >= seq_len {
            return None;
        }

        let actual_end = end.min(seq_len);
        let len = actual_end - start;

        if len < self.minimum_segment_length {
            return None;
        }

        let sequence: Vec<u8> = self
            .read
            .sequence()
            .as_ref()
            .iter()
            .skip(start)
            .take(len)
            .copied()
            .collect();

        let quality: Vec<u8> = self
            .read
            .quality_scores()
            .as_ref()
            .iter()
            .skip(start)
            .take(len)
            .map(|&q| q + 33)
            .collect();

        let metadata = SegmentMetadata::from_parts(
            self.name().ok()?.as_str(),
            self.viewpoint,
            ViewpointPosition::from_segment_type(segment_type),
            self.read_number(),
            self.flashed_status,
        );

        Some(Segment::<fastq::Record>::from_metadata(
            metadata, &sequence, &quality,
        ))
    }

    fn segments(&self) -> Option<Result<Vec<Segment<fastq::Record>>>> {
        if !self.is_viewpoint_read() {
            return None;
        }

        let positions = self.parse_segment_positions()?;
        let mut segments = Vec::new();

        for (seg_type, range) in positions {
            if let Some(segment) = self.extract_segment(seg_type, range) {
                segments.push(segment);
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
            flashed_status: FlashedStatus::Flashed,
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
            .build_from_path(self.bam_path.clone())
            .expect("Failed to build indexed reader");
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
            let query = reader
                .query(&header, &region)
                .expect("Failed to query region");
            for result in query.records() {
                let record = result?;
                counter += 1;
                if counter % 100_000 == 0 {
                    info!("Processed {} reads", counter);
                }

                let record_buf = noodles::sam::alignment::RecordBuf::try_from_alignment_record(
                    &header, &record,
                )?;
                let viewpoint_read = ViewpointRead::new(
                    region.name().to_str().unwrap(),
                    record_buf,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_segment_positions() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::{Cigar, Sequence};
        use noodles::sam::alignment::RecordBuf;

        // Construct a read with: 10S (Left) 50M (Viewpoint) 20S (Right)
        let cigar = Cigar::from(vec![
            Op::new(Kind::SoftClip, 10),
            Op::new(Kind::Match, 50),
            Op::new(Kind::SoftClip, 20),
        ]);

        let read = RecordBuf::builder()
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 80]))
            .build();

        let viewpoint_read = ViewpointRead {
            viewpoint: "test-vp",
            read,
            flashed_status: FlashedStatus::Flashed,
            minimum_segment_length: 5,
        };

        let positions = viewpoint_read.parse_segment_positions().unwrap();

        assert_eq!(positions.left(), (0, 10));
        assert_eq!(positions.viewpoint(), (10, 60));
        assert_eq!(positions.right(), (60, 80));
    }
}
