use anyhow::anyhow;
use anyhow::{Context, Result};
use bio::bio_types::strand::Strand;
use bio::io::fasta::Sequence;
use bio::utils;
use bstr::ByteSlice;
use flate2;
use itertools::{self, Itertools};
use log::{info, warn};
use noodles::fastq;
use noodles::fastq::record::Definition;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::header::record::value::map::header;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::de;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use crate::utils::{FlashedStatus, SegmentMetadata, SegmentType, ViewpointPosition};

/// Returns the strand type based on the reverse complement flag.
fn get_strand(is_reverse: bool) -> bio::bio_types::strand::Strand {
    if is_reverse {
        bio::bio_types::strand::Strand::Reverse
    } else {
        bio::bio_types::strand::Strand::Forward
    }
}

/// Returns the reference sequence ID as a Result to avoid unwrap().
fn get_reference_id(read: &noodles::bam::Record) -> Result<usize> {
    let id = read
        .reference_sequence_id()
        .ok_or_else(|| anyhow!("Missing reference sequence ID"))??;
    Ok(id)
}

/// Determines ligation junction positions while ensuring no unwrap().
fn get_ligation_positions(
    reporter: &noodles::bam::Record,
    capture: &noodles::bam::Record,
    segment: SegmentType,
    reporter_strand: bio::bio_types::strand::Strand,
    capture_strand: bio::bio_types::strand::Strand,
) -> Result<(usize, usize)> {
    let reporter_start = reporter
        .alignment_start()
        .ok_or_else(|| anyhow!("Missing reporter alignment start"))??
        .get();

    let capture_start = capture
        .alignment_start()
        .ok_or_else(|| anyhow!("Missing capture alignment start"))??
        .get();

    let reporter_end = reporter_start + reporter.sequence().len();
    let capture_end = capture_start + capture.sequence().len();

    match (segment, reporter_strand, capture_strand) {
        (SegmentType::LEFT, Strand::Forward, Strand::Forward) => Ok((reporter_end, capture_start)),
        (SegmentType::LEFT, Strand::Reverse, Strand::Reverse) => Ok((reporter_start, capture_end)),
        (SegmentType::LEFT, Strand::Forward, Strand::Reverse) => Ok((reporter_end, capture_end)),
        (SegmentType::LEFT, Strand::Reverse, Strand::Forward) => {
            Ok((reporter_start, capture_start))
        }
        (SegmentType::RIGHT, Strand::Forward, Strand::Forward) => Ok((reporter_start, capture_end)),
        (SegmentType::RIGHT, Strand::Reverse, Strand::Reverse) => Ok((reporter_end, capture_start)),
        (SegmentType::RIGHT, Strand::Forward, Strand::Reverse) => {
            Ok((reporter_start, capture_start))
        }
        (SegmentType::RIGHT, Strand::Reverse, Strand::Forward) => Ok((reporter_end, capture_end)),
        _ => Err(anyhow!(
            "Could not determine ligation junctions for given strands"
        )),
    }
}

pub struct MCCReadGroup {
    reads: Vec<noodles::bam::record::Record>,
    flashed_status: FlashedStatus,
}

pub struct PairsRecord {
    viewpoint_id: String,
    read_id: String,
    chr1: usize,
    pos1: usize,
    chr2: usize,
    pos2: usize,
    strand1: String,
    strand2: String,
}

impl PairsRecord {
    pub fn new(
        viewpoint_id: String,
        read_id: String,
        chr1: usize,
        pos1: usize,
        chr2: usize,
        pos2: usize,
        strand1: String,
        strand2: String,
    ) -> Self {
        // // Check that chromosome 1 occurs before chromosome 2 if not swap them
        // let (chr1, pos1, strand1, chr2, pos2, strand2) = if chr1 > chr2 {
        //     (chr2, pos2, strand2, chr1, pos1, strand1)
        // } else {
        //     (chr1, pos1, strand1, chr2, pos2, strand2)
        // };

        // // Check that pos1 is less than pos2 if not swap them
        // let (pos1, strand1, pos2, strand2) = if pos1 > pos2 {
        //     (pos2, strand2, pos1, strand1)
        // } else {
        //     (pos1, strand1, pos2, strand2)
        // };

        PairsRecord {
            viewpoint_id,
            read_id,
            chr1,
            pos1,
            chr2,
            pos2,
            strand1,
            strand2,
        }
    }
}

impl MCCReadGroup {
    pub fn new(reads: Vec<noodles::bam::record::Record>, flashed_status: FlashedStatus) -> Self {
        MCCReadGroup {
            reads,
            flashed_status,
        }
    }

    pub fn viewpoint_reads(&self) -> impl Iterator<Item = &noodles::bam::Record> {
        self.reads.iter().filter(|read| {
            let read_name = match read.name() {
                Some(name) => name.to_str().unwrap(),
                None => return false,
            };

            let read_name = SegmentMetadata::new(read_name);
            read_name.viewpoint_position() != ViewpointPosition::ALL
        })
    }

    pub fn contains_viewpoint(&self) -> bool {
        self.viewpoint_reads().count() > 0
    }

    pub fn any_mapped(&self) -> bool {
        self.reads.iter().any(|read| !read.flags().is_unmapped())
    }

    pub fn mapped_reads(&self) -> Vec<&noodles::bam::Record> {
        let reads = self
            .reads
            .iter()
            .filter(|read| !read.flags().is_unmapped())
            .collect();
        reads
    }

    pub fn reporters(&self) -> Vec<&noodles::bam::Record> {
        let has_viewpoint_read = self.contains_viewpoint();
        let mut reads = Vec::new();

        for read in &self.reads {
            let name = SegmentMetadata::from_read_name(read.name());
            let is_mapped = !read.flags().is_unmapped();
            let is_viewpoint = match name.viewpoint_position() {
                ViewpointPosition::ALL => true,
                _ => false,
            };

            if is_mapped && !is_viewpoint && has_viewpoint_read {
                reads.push(read);
            }
        }

        reads
    }

    pub fn captures(&self) -> Vec<&noodles::bam::Record> {
        let mut viewpoint_reads = self.viewpoint_reads().collect::<Vec<_>>();

        if viewpoint_reads.len() > 1 && self.flashed_status == FlashedStatus::FLASHED {
            // If the viewpoint is flashed, we only expect one capture read per viewpoint read
            // If there are more than one, we need to filter out the one with the highest mapping quality
            viewpoint_reads.sort_by_key(|read| {
                let qual = match read.mapping_quality() {
                    Some(qual) => qual.get() as i8,
                    None => 0,
                };

                qual * -1
            });
            viewpoint_reads.truncate(1);
        }

        viewpoint_reads
    }

    pub fn filter_mapped(&self) -> MCCReadGroup {
        MCCReadGroup::new(
            self.mapped_reads().into_iter().cloned().collect(),
            self.flashed_status,
        )
    }

    fn ligation_junctions(&self) -> Result<Vec<PairsRecord>> {
        let reporters = self.reporters();
        let captures = self.captures();
        let capture = captures
            .get(0)
            .ok_or_else(|| anyhow!("No capture read found"))?;

        let mut pairs = Vec::new();

        for reporter in reporters {
            let reporter_meta = SegmentMetadata::from_read_name(reporter.name());
            let reporter_segment =
                SegmentType::from_viewpoint_position(reporter_meta.viewpoint_position());
            let reporter_strand = get_strand(reporter.flags().is_reverse_complemented());
            let capture_strand = get_strand(capture.flags().is_reverse_complemented());

            let (pos1, pos2) = get_ligation_positions(
                &reporter,
                &capture,
                reporter_segment,
                reporter_strand,
                capture_strand,
            )?;

            let pairs_record = PairsRecord::new(
                reporter_meta.viewpoint().to_string(),
                reporter_meta.to_string(),
                get_reference_id(&reporter)?,
                pos1,
                get_reference_id(capture)?,
                pos2,
                reporter_strand.to_string(),
                capture_strand.to_string(),
            );

            pairs.push(pairs_record);
        }

        Ok(pairs)
    }
}

impl std::fmt::Display for MCCReadGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ReadGroup(\n{}\n)",
            self.reads
                .iter()
                .map(|read| format!("{:?}", read))
                .collect::<Vec<_>>()
                .join("\n")
        )
    }
}

impl std::fmt::Debug for MCCReadGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt(f)
    }
}

pub fn split_reads(bam: &str, output_directory: &str) -> Result<()> {
    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;
    let header = bam.read_header()?;
    let mut handles = HashMap::new();

    let mcc_groups = bam.records().into_iter().chunk_by(|r| match r {
        Ok(record) => SegmentMetadata::from_read_name(record.name())
            .parent_id()
            .to_string(),
        Err(_) => "UNKNOWN".to_string(),
    });

    for (_, reads) in mcc_groups.into_iter() {
        let reads = reads.collect::<Result<Vec<_>, _>>()?;
        let read_group = MCCReadGroup::new(reads, FlashedStatus::FLASHED);

        if read_group.contains_viewpoint() && read_group.any_mapped() {
            let read_group = read_group.filter_mapped();

            for reporter in read_group.reporters() {
                let viewpoint = SegmentMetadata::from_read_name(reporter.name())
                    .viewpoint()
                    .to_string();

                let handle = handles.entry(viewpoint.clone()).or_insert_with(|| {
                    let path = Path::new(output_directory).join(format!("{}.bam", &viewpoint));
                    let mut writer = noodles::bam::io::writer::Builder::default()
                        .build_from_path(path)
                        .expect("Could not create BAM writer");
                    writer
                        .write_header(&header)
                        .expect("Could not write header");
                    writer
                });

                handle.write_record(&header, reporter)?;
            }
        }
    }

    Ok(())
}

pub fn identify_ligation_junctions(bam: &str, output_directory: &str) -> Result<()> {
    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;
    let header = bam.read_header()?;
    let mut handles = HashMap::new();

    let ref_id_to_chromosome = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(ii, (chrom_name, chrom_map))| {
            let chrom_name = chrom_name.to_string();
            (ii, chrom_name)
        })
        .collect::<HashMap<_, _>>();

    let mcc_groups = bam.records().into_iter().chunk_by(|r| match r {
        Ok(record) => SegmentMetadata::from_read_name(record.name())
            .parent_id()
            .to_string(),
        Err(_) => "UNKNOWN".to_string(),
    });

    for (_, reads) in mcc_groups.into_iter() {
        let reads = reads.collect::<Result<Vec<_>, _>>()?;
        let read_group = MCCReadGroup::new(reads, FlashedStatus::FLASHED);

        if read_group.contains_viewpoint() && read_group.any_mapped() {
            let read_group = read_group.filter_mapped();
            let pairs = read_group.ligation_junctions()?;

            for pair in pairs {
                let handle = handles.entry(pair.viewpoint_id.clone()).or_insert_with(|| {
                    let path =
                        Path::new(output_directory).join(format!("{}.pairs", &pair.viewpoint_id));
                    let file = std::fs::File::create(path).expect("Could not create file");
                    let writer = std::io::BufWriter::new(file);
                    writer
                });

                let chrom_1 = ref_id_to_chromosome
                    .get(&pair.chr1)
                    .ok_or_else(|| anyhow!("Failed to get chromosome name"))?;
                let chrom_2 = ref_id_to_chromosome
                    .get(&pair.chr2)
                    .ok_or_else(|| anyhow!("Failed to get chromosome name"))?;

                writeln!(
                    handle,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    pair.read_id,
                    chrom_1,
                    pair.pos1,
                    chrom_2,
                    pair.pos2,
                    pair.strand1,
                    pair.strand2
                )
                .context("Could not write record")?;
            }
        }
    }

    Ok(())
}
