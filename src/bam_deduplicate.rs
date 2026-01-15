use std::collections::HashSet;
use std::path::Path;
use anyhow::{Context, Result, anyhow};
use noodles::sam;
use itertools::Itertools;
use pyo3::prelude::*;
use log::info;

use crate::utils::{SegmentMetadata, FlashedStatus};
use crate::mcc_data_handler::MCCReadGroup;
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

#[derive(Eq, PartialEq, Hash, Debug)]
struct MoleculeKey {
    coords: Vec<SegmentCoord>,
}

impl MoleculeKey {
    fn from_read_group(group: &MCCReadGroup) -> Self {
        let mut coords = Vec::new();
        for read in &group.reads {
            if let (Some(rid), Some(start_pos)) = (read.reference_sequence_id(), read.alignment_start()) {
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
        // Sort coordinates to make the key invariant to read order in the group
        coords.sort();
        Self { coords }
    }
}

pub fn deduplicate_bam(bam_path: &str, out_path: &str) -> Result<BamDeduplicationStats> {
    let mut reader = noodles::bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let mut writer = noodles::bam::io::writer::Builder::default().build_from_path(out_path)?;
    writer.write_header(&header)?;

    let mut stats = BamDeduplicationStats::new();
    let mut seen_molecules = HashSet::new();

    let mcc_groups = reader.records().into_iter().chunk_by(|r| {
        r.as_ref()
            .map(|record| SegmentMetadata::from_read_name(record.name()).parent_id().to_string())
            .unwrap_or_else(|_| "UNKNOWN".to_string())
    });

    for (parent_id, reads) in mcc_groups.into_iter() {
        if parent_id == "UNKNOWN" {
            // How to handle unknown groups? Usually skip or pass through.
            // For now, let's just skip them to be safe or pass them through.
            // Let's pass them through but they won't be deduplicated correctly.
            continue;
        }

        let reads_raw = reads.collect::<Result<Vec<_>, _>>()?;
        let mut reads_buf = Vec::new();
        for r in reads_raw {
            reads_buf.push(noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &r)?);
        }
        let group = MCCReadGroup::new(reads_buf, FlashedStatus::FLASHED);
        
        stats.total_molecules += 1;
        let key = MoleculeKey::from_read_group(&group);

        if seen_molecules.insert(key) {
            stats.unique_molecules += 1;
            // Write all reads in the group
            for read in group.reads {
                writer.write_alignment_record(&header, &read)?;
            }
        } else {
            stats.duplicate_molecules += 1;
        }
    }

    writer.try_finish()?;
    info!("BAM deduplication complete: {} total, {} unique, {} duplicates", 
          stats.total_molecules, stats.unique_molecules, stats.duplicate_molecules);

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::bam;

    #[test]
    fn test_segment_coord_sorting() {
        let c1 = SegmentCoord { chrom_id: 1, start: 100, end: 200, is_reverse: false };
        let c2 = SegmentCoord { chrom_id: 1, start: 50, end: 150, is_reverse: false };
        let c3 = SegmentCoord { chrom_id: 0, start: 100, end: 200, is_reverse: false };
        
        let mut coords = vec![c1, c2, c3];
        coords.sort();
        
        assert_eq!(coords[0].chrom_id, 0);
        assert_eq!(coords[1].start, 50);
        assert_eq!(coords[2].start, 100);
    }

    #[test]
    fn test_molecule_key_invariance() {
        use noodles::sam::alignment::record_buf::{Sequence, QualityScores};
        use noodles::sam::alignment::RecordBuf;
        use noodles::core::Position;

        // MoleculeKey should be invariant to the order of reads in the group
        let r1 = RecordBuf::builder()
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::try_from(100).unwrap())
            .set_sequence(Sequence::from(vec![b'A'; 50]))
            .set_quality_scores(QualityScores::from(vec![b'!'; 50]))
            .build();

        let r2 = RecordBuf::builder()
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::try_from(200).unwrap())
            .set_sequence(Sequence::from(vec![b'C'; 50]))
            .set_quality_scores(QualityScores::from(vec![b'!'; 50]))
            .build();

        let group1 = MCCReadGroup::new(vec![r1.clone(), r2.clone()], FlashedStatus::FLASHED);
        let group2 = MCCReadGroup::new(vec![r2, r1], FlashedStatus::FLASHED);

        let key1 = MoleculeKey::from_read_group(&group1);
        let key2 = MoleculeKey::from_read_group(&group2);

        assert_eq!(key1, key2);
    }
}
