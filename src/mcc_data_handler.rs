use crate::utils::{FlashedStatus, SegmentMetadata, SegmentType, ViewpointPosition};
use anyhow::{anyhow, Context, Result};
use bio::bio_types::strand::Strand;
use bio::io::fasta::Sequence;
use bio::utils;
use bstr::ByteSlice;
use flate2;
use itertools::{self, Itertools};
use log::{info, warn};
use noodles::fastq;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::header::record::value::map::header;
use noodles::sam::header::record::value::{map::ReadGroup, Map};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::de;
use std::collections::{HashMap, HashSet};
use std::default;
use std::io::BufRead;
use std::io::Write; // Import the Write trait for writeln! macro
use std::path::{Path, PathBuf};

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
                reporter_meta.viewpoint_name().to_string(),
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


    pub fn is_valid(&self, chrom1_length: usize, chrom2_length:usize) -> bool {
        // Check that the positions are within the chromosome lengths
        if self.pos1 > chrom1_length || self.pos2 > chrom2_length {
            return false;
        } else if self.pos1 == 0 || self.pos2 == 0 {
            return false;
        }
        true
    }


}

pub fn annotate_bam(bam: &str, out: &str) -> Result<()> {
    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;
    let header = bam.read_header()?;

    let temp = PathBuf::from(out).with_extension("temp.bam");
    if temp.exists() {
        std::fs::remove_file(&temp).context("Could not remove existing file")?;
    }
    let mut writer = noodles::bam::io::writer::Builder::default().build_from_path(&temp)?;
    writer.write_header(&header)?;

    // Store read groups in a hash set
    let mut read_groups_set = HashSet::new();

    // Viewpoint tag
    let vp_tag = noodles::sam::alignment::record::data::field::Tag::new(b'V', b'P');

    // OC tag -- oligo coordinate tag
    let oc_tag = noodles::sam::alignment::record::data::field::Tag::new(b'O', b'C');

    // Reporter tag
    let reporter_tag = noodles::sam::alignment::record::data::field::Tag::new(b'R', b'T');

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

            // Get the details of the viewpoint
            let example_read = read_group.reads.get(0).context("No reads in read group")?;
            let viewpoint = SegmentMetadata::from_read_name(example_read.name())
                .viewpoint()
                .to_string();

            // Get the viewpoint name
            let viewpoint_name = viewpoint
                .split_once("-")
                .context("Could not split viewpoint name")?
                .0;

            // Get the oligo coordinate from the read name
            let oligo_coordinate = SegmentMetadata::from_read_name(example_read.name())
                .oligo_coordinates()
                .to_string();

            read_groups_set.insert(viewpoint_name.to_string());

            // Write the capture reads to the output
            for capture_read in read_group.captures() {
                let mut record_sam = noodles::sam::alignment::RecordBuf::try_from_alignment_record(
                    &header,
                    capture_read,
                )?;

                // Add the reporter tag to the record
                record_sam.data_mut().insert(reporter_tag, Value::Int8(0));

                // Add the OC tag to the record
                record_sam
                    .data_mut()
                    .insert(oc_tag, Value::String(oligo_coordinate.clone().into()));
                // Add the VP tag to the record -- this is just the viewpoint name before the first "-"
                record_sam
                    .data_mut()
                    .insert(vp_tag, Value::String(viewpoint_name.into()));
                // Add the read group tag to the record
                record_sam.data_mut().insert(
                    noodles::sam::alignment::record::data::field::Tag::READ_GROUP,
                    Value::String(viewpoint_name.into()),
                );

                noodles::sam::alignment::io::Write::write_alignment_record(
                    &mut writer,
                    &header,
                    &record_sam,
                )?;
            }

            // Write the reporter reads to the output
            for reporter in read_group.reporters() {
                let mut record_sam = noodles::sam::alignment::RecordBuf::try_from_alignment_record(
                    &header, reporter,
                )?;

                // Add the reporter tag to the record
                record_sam.data_mut().insert(reporter_tag, Value::Int8(1));

                // Add the OC tag to the record
                record_sam
                    .data_mut()
                    .insert(oc_tag, Value::String(oligo_coordinate.clone().into()));

                // Add the VP tag to the record -- this is just the viewpoint name before the first "-"
                record_sam
                    .data_mut()
                    .insert(vp_tag, Value::String(viewpoint_name.into()));

                // Add the read group tag to the record
                record_sam.data_mut().insert(
                    noodles::sam::alignment::record::data::field::Tag::READ_GROUP,
                    Value::String(viewpoint_name.into()),
                );

                noodles::sam::alignment::io::Write::write_alignment_record(
                    &mut writer,
                    &header,
                    &record_sam,
                )?;
            }
        }
    }

    // Close the writer
    writer.try_finish()?;

    // Need to re-header the output file to include the read groups
    let mut bam_in = noodles::bam::io::reader::Builder::default().build_from_path(&temp)?;
    let mut header = bam_in.read_header()?;

    for rg in read_groups_set {
        let read_group = Map::<ReadGroup>::default();
        header.read_groups_mut().insert(rg.into(), read_group);
    }
    let mut bam_out = noodles::bam::io::writer::Builder::default()
        .build_from_path(out)
        .context("Could not create output file")?;

    bam_out.write_header(&header).context("Could not write header")?;

    std::io::copy(bam_in.get_mut(), bam_out.get_mut())?;

    // Remove the temporary file
    std::fs::remove_file(&temp).context("Could not remove temporary file")?;
    info!("Finished annotating BAM file");
    info!("Output written to {}", out);

    Ok(())
}

pub fn identify_ligation_junctions(bam: &str, output_directory: &str) -> Result<()> {
    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;
    let header = bam.read_header()?;
    let mut handles = HashMap::new();

    // Will ideally depend on the bamnado package here. This is not published yet

    // Reference ID to chromosome name mapping
    // This is a temporary solution until we can use the bamnado BamStats structure
    let ref_id_to_chromosome = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(ii, (chrom_name, chrom_map))| {
            let chrom_name = chrom_name.to_string();
            (ii, chrom_name)
        })
        .collect::<HashMap<_, _>>();

    // Get Chromosome lengths
    let chrom_lengths = header
        .reference_sequences()
        .iter()
       .map(|(chrom_name, chrom_map)| {
            let chrom_name = chrom_name.to_string();
            let chrom_length = chrom_map.length().into();
            (chrom_name, chrom_length)
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

                // Check if the pair is valid
                let chrom1_length = chrom_lengths
                    .get(chrom_1)
                    .ok_or_else(|| anyhow!("Chromosome 1 not found"))?;
                let chrom2_length = chrom_lengths
                    .get(chrom_2)
                    .ok_or_else(|| anyhow!("Chromosome 2 not found"))?;

                if !pair.is_valid(*chrom1_length, *chrom2_length) {
                    continue;
                }

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
