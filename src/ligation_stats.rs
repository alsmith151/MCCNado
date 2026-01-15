use anyhow::{Context, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{BufWriter, Write};

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

struct LigationTags {
    viewpoint: Tag,
    oligo_coordinate: Tag,
    reporter: Tag,
}

impl LigationTags {
    fn new() -> Self {
        Self {
            viewpoint: Tag::new(b'V', b'P'),
            oligo_coordinate: Tag::new(b'O', b'C'),
            reporter: Tag::new(b'R', b'T'),
        }
    }
}

fn extract_tag_string(
    record_data: &noodles::sam::alignment::record_buf::Data,
    tag: Tag,
) -> Option<String> {
    match record_data.get(&tag) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

fn is_reporter_read(record_data: &noodles::sam::alignment::record_buf::Data, tag: Tag) -> bool {
    match record_data.get(&tag) {
        Some(Value::UInt8(v)) => *v != 0,
        Some(Value::Int8(v)) => *v != 0,
        _ => false,
    }
}

fn write_stats_json(output_path: &str, stats: &HashMap<String, LigationStats>) -> Result<()> {
    let file = std::fs::File::create(output_path).context("Failed to create output file")?;
    let mut writer = BufWriter::new(file);
    let json = serde_json::to_string(stats).context("Failed to serialize ligation stats")?;
    writer
        .write_all(json.as_bytes())
        .context("Failed to write ligation stats to file")?;
    Ok(())
}

pub fn get_ligation_stats(bam_path: &str, output_path: &str) -> Result<()> {
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let mut ligation_stats: HashMap<String, LigationStats> = HashMap::new();

    let refid_to_name: HashMap<usize, String> = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (i, name.to_string()))
        .collect();

    let tags = LigationTags::new();

    for record in reader.records() {
        let record = noodles::sam::alignment::record_buf::RecordBuf::try_from_alignment_record(
            &header, &record?,
        )?;
        let data = record.data();

        if !is_reporter_read(data, tags.reporter) {
            continue;
        }

        let viewpoint_name = extract_tag_string(data, tags.viewpoint).context("Missing VP tag")?;
        let oligo_coordinate =
            extract_tag_string(data, tags.oligo_coordinate).context("Missing OC tag")?;

        let entry = ligation_stats
            .entry(viewpoint_name)
            .or_insert_with(LigationStats::new);
        entry.n_total += 1;

        let chrom_read_id = record
            .reference_sequence_id()
            .context("Missing reference sequence ID")?;
        let chrom_read = refid_to_name
            .get(&chrom_read_id)
            .context("Could not get chromosome name")?;

        let chrom_viewpoint = oligo_coordinate
            .split_once('-')
            .context("Could not split oligo coordinate")?
            .0;

        if chrom_read == chrom_viewpoint {
            entry.n_cis += 1;
        } else {
            entry.n_trans += 1;
        }
    }

    write_stats_json(output_path, &ligation_stats)?;
    Ok(())
}
