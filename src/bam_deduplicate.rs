use anyhow::Result;
use itertools::Itertools;
use log::info;
use noodles::sam::alignment::record::data::field::Tag;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hasher;
use twox_hash::XxHash64;

use crate::utils::SegmentMetadata;
use noodles::sam::alignment::io::Write;

#[derive(Debug, serde::Serialize, serde::Deserialize, Clone)]
#[pyclass]
pub struct BamDeduplicationStats {
    #[pyo3(get)]
    pub total_molecules: u64,
    #[pyo3(get)]
    pub unique_molecules: u64,
    #[pyo3(get)]
    pub duplicate_molecules: u64,
}

impl BamDeduplicationStats {
    fn new() -> Self {
        Self {
            total_molecules: 0,
            unique_molecules: 0,
            duplicate_molecules: 0,
        }
    }
}

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Debug)]
struct SegmentCoord {
    chrom_id: usize,
    start: usize,
    end: usize,
    is_reverse: bool,
}

// Compute key directly from raw records without intermediate allocations
fn compute_molecule_key_from_raw(reads: &[noodles::bam::Record]) -> Result<u64> {
    let mut coords: Vec<SegmentCoord> = Vec::with_capacity(reads.len());

    for read in reads {
        if let (Some(rid), Some(start_pos)) = (read.reference_sequence_id(), read.alignment_start())
        {
            let rid = rid?; // Handle Result
            let start_pos = start_pos?;
            let start = start_pos.get();
            let end = start + read.sequence().len();
            coords.push(SegmentCoord {
                chrom_id: rid,
                start,
                end,
                is_reverse: read.flags().is_reverse_complemented(),
            });
        }
    }

    coords.sort();
    coords.dedup(); // Fix: Deduplicate segments to handle RT=0 vs RT=1 redundancy

    let mut hasher = XxHash64::with_seed(0);
    for coord in &coords {
        hasher.write_usize(coord.chrom_id);
        hasher.write_usize(coord.start);
        hasher.write_usize(coord.end);
        hasher.write_u8(coord.is_reverse as u8);
    }
    Ok(hasher.finish())
}

fn filter_output_reads(
    reads: Vec<noodles::sam::alignment::RecordBuf>,
) -> Vec<noodles::sam::alignment::RecordBuf> {
    // Map of ReadName -> (Index in original list, Record)
    // We want to keep RT=1 over RT=0 if names match
    let mut unique_reads: HashMap<String, (usize, noodles::sam::alignment::RecordBuf)> =
        HashMap::with_capacity(reads.len());

    for (i, read) in reads.into_iter().enumerate() {
        let name = read.name().map(|n| n.to_string()).unwrap_or_default();

        // Check for RT tag: 'RT' = 'R' 'T'
        let is_reporter = read
            .data()
            .get(&Tag::new(b'R', b'T'))
            .and_then(|v| v.as_int())
            .map(|v| v == 1)
            .unwrap_or(false);

        match unique_reads.get_mut(&name) {
            Some((_, existing)) => {
                // If we found a duplicate, check if new one is reporter and old one isn't
                let existing_is_reporter = existing
                    .data()
                    .get(&Tag::new(b'R', b'T'))
                    .and_then(|v| v.as_int())
                    .map(|v| v == 1)
                    .unwrap_or(false);

                if is_reporter && !existing_is_reporter {
                    *existing = read;
                }
            }
            None => {
                unique_reads.insert(name, (i, read));
            }
        }
    }

    // Restore original order
    let mut result: Vec<_> = unique_reads.into_values().collect();
    result.sort_by_key(|(i, _)| *i);
    result.into_iter().map(|(_, r)| r).collect()
}

pub fn deduplicate_bam(bam_path: &str, out_path: &str) -> Result<BamDeduplicationStats> {
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let mut writer = noodles::bam::io::writer::Builder.build_from_path(out_path)?;
    writer.write_header(&header)?;

    let mut stats = BamDeduplicationStats::new();
    // Pre-size with a reasonable guess to avoid early reallocations
    let mut seen_molecules = HashSet::with_capacity(100_000);

    let mcc_groups = reader.records().chunk_by(|r| {
        r.as_ref()
            .map(|record| {
                SegmentMetadata::from_read_name(record.name())
                    .parent_id()
                    .to_string()
            })
            .unwrap_or_else(|_| "UNKNOWN".to_string())
    });

    for (parent_id, reads) in mcc_groups.into_iter() {
        if parent_id == "UNKNOWN" {
            // How to handle unknown groups? Usually skip or pass through.
            // For now, let's just skip them to be safe or pass them through.
            // Let's pass them through but they won't be deduplicated correctly.
            continue;
        }

        let reads_raw: Vec<_> = reads.collect::<Result<Vec<_>, _>>()?;
        stats.total_molecules += 1;

        // Compute key from raw records (no RecordBuf allocation needed)
        let key = compute_molecule_key_from_raw(&reads_raw)?;

        if seen_molecules.insert(key) {
            stats.unique_molecules += 1;

            // Only now convert to RecordBuf for writing
            let reads_buf: Vec<_> = reads_raw
                .iter()
                .map(|r| noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, r))
                .collect::<Result<_, _>>()?;

            // Filter duplicates from output (removing RT=0/RT=1 redundancy)
            let filtered_reads = filter_output_reads(reads_buf);

            for read in filtered_reads {
                writer.write_alignment_record(&header, &read)?;
            }
        } else {
            stats.duplicate_molecules += 1;
        }
    }

    writer.try_finish()?;
    info!(
        "BAM deduplication complete: {} total, {} unique, {} duplicates",
        stats.total_molecules, stats.unique_molecules, stats.duplicate_molecules
    );

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_coord_sorting() {
        let c1 = SegmentCoord {
            chrom_id: 1,
            start: 100,
            end: 200,
            is_reverse: false,
        };
        let c2 = SegmentCoord {
            chrom_id: 1,
            start: 50,
            end: 150,
            is_reverse: false,
        };
        let c3 = SegmentCoord {
            chrom_id: 0,
            start: 100,
            end: 200,
            is_reverse: false,
        };

        let mut coords = [c1, c2, c3];
        coords.sort();

        assert_eq!(coords[0].chrom_id, 0);
        assert_eq!(coords[1].start, 50);
        assert_eq!(coords[2].start, 100);
    }
}
