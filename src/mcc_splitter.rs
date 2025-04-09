use anyhow::{Context, Result};
use bio::io::fasta::Sequence;
use bio::utils;
use bstr::ByteSlice;
use flate2;
use itertools::{self, Itertools};
use log::{info, warn};
use noodles::fastq;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::header::record::value::map::header;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::de;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead};
use std::path::{Path, PathBuf};

use crate::utils::{FlashedStatus, SegmentMetadata, ViewpointPosition};

pub struct MCCReadGroup {
    reads: Vec<noodles::bam::record::Record>,
    flashed_status: FlashedStatus,
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
                    writer.write_header(&header).expect("Could not write header");
                    writer
                });

                handle.write_record(&header, reporter)?;
            }
        }
    }

    Ok(())
}


pub fn add_viewpoint_tag(
    bam: &str,
) -> Result<()> {


    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;
    let header = bam.read_header()?;

    let stdout = std::io::stdout().lock();
    let mut writer = noodles::sam::io::Writer::new(stdout);
    writer.write_header(&header)?;

    // Viewpoint tag
    let vp_tag = noodles::sam::alignment::record::data::field::Tag::new(b'V', b'P');

    // OC tag -- oligo coordinate tag
    let oc_tag = noodles::sam::alignment::record::data::field::Tag::new(b'O', b'C');

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

                let mut record_sam = noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, reporter)?;

                // Add the OC tag to the record
                record_sam.data_mut().insert(oc_tag, Value::String(viewpoint.clone().into()));
                
                // Add the VP tag to the record -- this is just the viewpoint name before the first "-"
                let vp_name = viewpoint.split_once("-").context("Could not split viewpoint name")?.0;
                record_sam.data_mut().insert(vp_tag, Value::String(vp_name.into()));
                

                // Add the read group tag to the record
                record_sam.data_mut().insert(
                    noodles::sam::alignment::record::data::field::Tag::READ_GROUP,
                    Value::String(
                        vp_name.into()
                    )
                );

                writer.write_alignment_record(&header, &record_sam)?;


            }
        }
    }







    Ok(())







}