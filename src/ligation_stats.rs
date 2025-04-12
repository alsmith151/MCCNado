use anyhow::{Context, Result};
use bio::io::fasta::Sequence;
use bio::utils;
use bstr::ByteSlice;
use flate2;
use itertools::{self, Itertools};
use log::{info, warn};
use noodles::fastq;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::header::record::value::map::header;
use pyo3::types::PyDict;
use serde::{de, Deserialize, Serialize};
use std::any::Any;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize, Clone, Deserialize)]
struct LigationStats {
    n_cis: u64,
    n_trans: u64,
    n_total: u64,
}

impl LigationStats {
    pub fn new() -> Self {
        Self {
            n_cis: 0,
            n_trans: 0,
            n_total: 0,
        }
    }
}

pub fn get_ligation_stats(bam: &str, output: &str) -> Result<()> {
    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;
    let header = bam.read_header()?;
    let mut ligation_stats: HashMap<String, LigationStats> = HashMap::new();

    let refid_to_chrom_name: HashMap<usize, String> = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, map))| {
            let name = name.to_string();
            (i, name)
        })
        .collect();

    // Viewpoint tag
    let vp_tag = noodles::sam::alignment::record::data::field::Tag::new(b'V', b'P');

    // OC tag -- oligo coordinate tag
    let oc_tag = noodles::sam::alignment::record::data::field::Tag::new(b'O', b'C');

    // Reporter tag
    let reporter_tag = noodles::sam::alignment::record::data::field::Tag::new(b'R', b'T');

    for record in bam.records() {
        let record = noodles::sam::alignment::record_buf::RecordBuf::try_from_alignment_record(
            &header, &record?,
        )?;
        let rdata = record.data();

        // Extract the viewpoint name from the read
        let viewpoint_name = rdata.get(&vp_tag).context("Missing VP tag")?;
        let viewpoint_name = match viewpoint_name {
            Value::String(s) => s.to_string(),
            _ => {
                warn!("Invalid VP tag value");
                continue;
            }
        };

        // Extract the oligo coordinate tag from the read
        let oligo_coordinate = rdata.get(&oc_tag).context("Missing OC tag")?;
        let oligo_coordinate = match oligo_coordinate {
            Value::String(s) => s.to_string(),
            _ => {
                warn!("Invalid OC tag value");
                continue;
            }
        };

        // Extract if the read is a repoter from the read
        let is_reporter = rdata.get(&reporter_tag).context("Missing RT tag")?;
        let is_reporter = match is_reporter {
            Value::UInt8(s) => {
                if *s == 0 {
                    false
                } else {
                    true
                }
            },
            Value::Int8(s) => {
                if *s == 0 {
                    false
                } else {
                    true
                }
            },
            _ => {
                warn!("Invalid RT tag value: {:?}", is_reporter);
                false
            }
        };

        if !is_reporter {
            continue;
        } else {
            let mut stat_entry = ligation_stats
                .entry(viewpoint_name)
                .or_insert(LigationStats::new());
            stat_entry.n_total += 1;

            let chromosome_read = record
                .reference_sequence_id()
                .context("Missing reference sequence ID")?;

            let chromosome_read = refid_to_chrom_name
                .get(&chromosome_read)
                .context("Could not get chromosome name")?
                .to_string();

            let chromosome_viewpoint = oligo_coordinate
                .split_once("-")
                .context("Could not split oligo coordinate")?
                .0;

            if chromosome_read == chromosome_viewpoint {
                stat_entry.n_cis += 1;
            } else {
                stat_entry.n_trans += 1;
            }
        }
    }

    // Write the ligation stats to the output file
    let stats_file = std::fs::File::create(output).context("Failed to create output file")?;

    let mut stats_writer = BufWriter::new(stats_file);

    let stats =
        serde_json::to_string(&ligation_stats).context("Failed to serialize ligation stats")?;

    stats_writer
        .write(stats.as_bytes())
        .context("Failed to write ligation stats to file")?;

    Ok(())
}
