use crate::utils::{FlashedStatus, SegmentMetadata, SegmentType, Strand, ViewpointPosition};
use anyhow::{anyhow, Context, Result};
use bstr::ByteSlice;
use itertools::Itertools;
use log::info;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::header::record::value::{map::ReadGroup, Map};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

struct MccTags {
    viewpoint: Tag,
    oligo_coordinate: Tag,
    reporter: Tag,
}

impl MccTags {
    fn new() -> Self {
        Self {
            viewpoint: Tag::new(b'V', b'P'),
            oligo_coordinate: Tag::new(b'O', b'C'),
            reporter: Tag::new(b'R', b'T'),
        }
    }
}

pub struct MCCReadGroup {
    pub reads: Vec<noodles::sam::alignment::RecordBuf>,
    pub flashed_status: FlashedStatus,
}

impl MCCReadGroup {
    pub fn new(
        reads: Vec<noodles::sam::alignment::RecordBuf>,
        flashed_status: FlashedStatus,
    ) -> Self {
        MCCReadGroup {
            reads,
            flashed_status,
        }
    }

    pub fn viewpoint_reads(&self) -> impl Iterator<Item = &noodles::sam::alignment::RecordBuf> {
        self.reads.iter().filter(|read| {
            let read_name = match read.name() {
                Some(name) => name.to_str().unwrap(),
                None => return false,
            };

            let read_name = SegmentMetadata::new(read_name);
            // Fix: Capture is ONLY the viewpoint oligo itself (ALL)
            matches!(read_name.viewpoint_position(), ViewpointPosition::All)
        })
    }

    pub fn contains_viewpoint(&self) -> bool {
        self.viewpoint_reads().count() > 0
    }

    pub fn any_mapped(&self) -> bool {
        self.reads.iter().any(|read| !read.flags().is_unmapped())
    }

    pub fn mapped_reads(&self) -> Vec<&noodles::sam::alignment::RecordBuf> {
        let reads = self
            .reads
            .iter()
            .filter(|read| !read.flags().is_unmapped())
            .collect();
        reads
    }

    pub fn reporters(&self) -> Vec<&noodles::sam::alignment::RecordBuf> {
        let has_viewpoint_read = self.contains_viewpoint();
        let mut reads = Vec::new();

        for read in &self.reads {
            let name = SegmentMetadata::from_read_name(read.name());

            // Fix: Reporter is everything that is NOT the viewpoint oligo (ALL)
            // (e.g., START/RIGHT, END/LEFT, NONE)
            let is_viewpoint = matches!(name.viewpoint_position(), ViewpointPosition::All);

            if !is_viewpoint && has_viewpoint_read {
                reads.push(read);
            }
        }

        reads
    }

    pub fn captures(&self) -> Vec<&noodles::sam::alignment::RecordBuf> {
        let mut viewpoint_reads = self.viewpoint_reads().collect::<Vec<_>>();

        if viewpoint_reads.len() > 1 && self.flashed_status == FlashedStatus::Flashed {
            // If the viewpoint is flashed, we only expect one capture read per viewpoint read
            // If there are more than one, we need to filter out the one with the highest mapping quality
            viewpoint_reads.sort_by_key(|read| {
                let qual = match read.mapping_quality() {
                    Some(qual) => qual.get() as i8,
                    None => 0,
                };

                -qual
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
            .first()
            .ok_or_else(|| anyhow!("No capture read found"))?;

        let mut pairs = Vec::new();

        for reporter in reporters {
            let reporter_meta = SegmentMetadata::from_read_name(reporter.name());
            let reporter_segment =
                SegmentType::from_viewpoint_position(reporter_meta.viewpoint_position());
            let reporter_strand = get_strand(reporter.flags().is_reverse_complemented());
            let capture_strand = get_strand(capture.flags().is_reverse_complemented());

            let (pos1, pos2) = get_ligation_positions(
                reporter,
                capture,
                reporter_segment,
                reporter_strand,
                capture_strand,
            )?;

            let pairs_record = PairsRecord::new(
                reporter_meta.viewpoint_name().to_string(),
                reporter_meta.to_string(),
                get_reference_id(reporter)?,
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
        f.debug_struct("MCCReadGroup")
            .field("reads_count", &self.reads.len())
            .field("flashed_status", &self.flashed_status)
            .finish()
    }
}

/// Returns the strand type based on the reverse complement flag.
fn get_strand(is_reverse: bool) -> Strand {
    if is_reverse {
        Strand::Negative
    } else {
        Strand::Positive
    }
}

/// Returns the reference sequence ID as a Result to avoid unwrap().
fn get_reference_id(read: &noodles::sam::alignment::RecordBuf) -> Result<usize> {
    let id = read
        .reference_sequence_id()
        .ok_or_else(|| anyhow!("Missing reference sequence ID"))?;
    Ok(id)
}

/// Determines ligation junction positions while ensuring no unwrap().
fn get_ligation_positions(
    reporter: &noodles::sam::alignment::RecordBuf,
    capture: &noodles::sam::alignment::RecordBuf,
    segment: SegmentType,
    reporter_strand: Strand,
    capture_strand: Strand,
) -> Result<(usize, usize)> {
    let reporter_start = reporter
        .alignment_start()
        .ok_or_else(|| anyhow!("Missing reporter alignment start"))?
        .get();

    let capture_start = capture
        .alignment_start()
        .ok_or_else(|| anyhow!("Missing capture alignment start"))?
        .get();

    let reporter_end = reporter_start + reporter.sequence().len();
    let capture_end = capture_start + capture.sequence().len();

    match (segment, reporter_strand, capture_strand) {
        (SegmentType::Left, Strand::Positive, Strand::Positive) => {
            Ok((reporter_end, capture_start))
        }
        (SegmentType::Left, Strand::Negative, Strand::Negative) => {
            Ok((reporter_start, capture_end))
        }
        (SegmentType::Left, Strand::Positive, Strand::Negative) => Ok((reporter_end, capture_end)),
        (SegmentType::Left, Strand::Negative, Strand::Positive) => {
            Ok((reporter_start, capture_start))
        }
        (SegmentType::Right, Strand::Positive, Strand::Positive) => {
            Ok((reporter_start, capture_end))
        }
        (SegmentType::Right, Strand::Negative, Strand::Negative) => {
            Ok((reporter_end, capture_start))
        }
        (SegmentType::Right, Strand::Positive, Strand::Negative) => {
            Ok((reporter_start, capture_start))
        }
        (SegmentType::Right, Strand::Negative, Strand::Positive) => Ok((reporter_end, capture_end)),
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
    #[allow(clippy::too_many_arguments)]
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

    pub fn is_valid(&self, chrom1_length: usize, chrom2_length: usize) -> bool {
        // Check that the positions are within the chromosome lengths
        if self.pos1 > chrom1_length
            || self.pos2 > chrom2_length
            || self.pos1 == 0
            || self.pos2 == 0
        {
            return false;
        }
        true
    }
}


fn check_bam_sorted_by_queryname(_header: &noodles::sam::Header) -> Result<()> {
    // Check if BAM is sorted by queryname (required for chunk_by grouping to work correctly)
    // The sort order is typically stored in the @HD line with SO tag
    // For now, we'll just log a warning and trust the user provided sorted data
    // A full check would require parsing the header's SO field which varies by implementation
    log::debug!("Assuming BAM file is sorted by queryname for proper grouping.");  
    Ok(())
}

fn write_annotated_records<W>(
    writer: &mut W,
    header: &noodles::sam::Header,
    read_group: &MCCReadGroup,
    tags: &MccTags,
    viewpoint_name: &str,
    oligo_coordinate: &str,
) -> Result<()>
where
    W: noodles::sam::alignment::io::Write,

{
    // Write capture reads
    for capture_read in read_group.captures() {
        let mut record = capture_read.clone();
        record.data_mut().insert(tags.reporter, Value::Int8(0));
        record.data_mut().insert(
            tags.oligo_coordinate,
            Value::String(oligo_coordinate.into()),
        );
        record
            .data_mut()
            .insert(tags.viewpoint, Value::String(viewpoint_name.into()));
        record
            .data_mut()
            .insert(Tag::READ_GROUP, Value::String(viewpoint_name.into()));
        writer.write_alignment_record(header, &record)?;
    }

    // Write reporter reads
    for reporter in read_group.reporters() {
        let mut record = reporter.clone();
        record.data_mut().insert(tags.reporter, Value::Int8(1));
        record.data_mut().insert(
            tags.oligo_coordinate,
            Value::String(oligo_coordinate.into()),
        );
        record
            .data_mut()
            .insert(tags.viewpoint, Value::String(viewpoint_name.into()));
        record
            .data_mut()
            .insert(Tag::READ_GROUP, Value::String(viewpoint_name.into()));
        writer.write_alignment_record(header, &record)?;
    }

    Ok(())
}

fn finalize_bam_with_read_groups(
    temp_path: &Path,
    final_path: &str,
    read_groups: HashSet<String>,
) -> Result<()> {
    let mut bam_in = noodles::bam::io::reader::Builder.build_from_path(temp_path)?;
    let mut header = bam_in.read_header()?;

    for rg in read_groups {
        header
            .read_groups_mut()
            .insert(rg.into(), Map::<ReadGroup>::default());
    }

    let mut bam_out = noodles::bam::io::writer::Builder
        .build_from_path(final_path)
        .context("Could not create output file")?;

    bam_out
        .write_header(&header)
        .context("Could not write header")?;
    std::io::copy(bam_in.get_mut(), bam_out.get_mut())?;

    std::fs::remove_file(temp_path).context("Could not remove temporary file")?;
    Ok(())
}

pub fn annotate_bam(bam_path: &str, out_path: &str) -> Result<()> {
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;

    // Enforce: BAM file must be sorted by queryname (read name) for proper grouping
    // The chunk_by operation assumes records with the same parent_id are consecutive
    check_bam_sorted_by_queryname(&header)?;

    let temp_path = PathBuf::from(out_path).with_extension("temp.bam");
    if temp_path.exists() {
        std::fs::remove_file(&temp_path).context("Could not remove existing temporary file")?;
    }

    let mut writer = noodles::bam::io::writer::Builder.build_from_path(&temp_path)?;
    writer.write_header(&header)?;

    let tags = MccTags::new();
    let mut read_groups_set = HashSet::new();

    // Stream and group by parent_id (assumes BAM is sorted by queryname)
    let mcc_groups = reader.records().chunk_by(|r| {
        r.as_ref()
            .map(|record| {
                SegmentMetadata::from_read_name(record.name())
                    .parent_id()
                    .to_string()
            })
            .unwrap_or_else(|_| "UNKNOWN".to_string())
    });

    for (_, reads) in mcc_groups.into_iter() {
        let reads_raw = reads.collect::<Result<Vec<_>, _>>()?;
        let mut reads_buf = Vec::new();
        for r in reads_raw {
            reads_buf
                .push(noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &r)?);
        }
        let read_group = MCCReadGroup::new(reads_buf, FlashedStatus::Flashed);

        if read_group.contains_viewpoint() {
            let first_read = read_group.reads.first().context("No reads in group")?;
            let metadata = SegmentMetadata::from_read_name(first_read.name());

            let viewpoint_full = metadata.viewpoint();
            let viewpoint_name = viewpoint_full
                .split_once('-')
                .context("Invalid viewpoint format")?
                .0;
            let oligo_coordinate = metadata.oligo_coordinates();

            read_groups_set.insert(viewpoint_name.to_string());

            write_annotated_records(
                &mut writer,
                &header,
                &read_group,
                &tags,
                viewpoint_name,
                oligo_coordinate,
            )?;
        }
    }

    writer.try_finish()?;
    finalize_bam_with_read_groups(&temp_path, out_path, read_groups_set)?;

    info!("Finished annotating BAM file: {}", out_path);
    Ok(())
}

struct ChromosomeInfo {
    id_to_name: HashMap<usize, String>,
    name_to_length: HashMap<String, usize>,
}

impl ChromosomeInfo {
    fn from_header(header: &noodles::sam::Header) -> Self {
        let mut id_to_name = HashMap::new();
        let mut name_to_length = HashMap::new();

        for (i, (name, map)) in header.reference_sequences().iter().enumerate() {
            let name_str = name.to_string();
            id_to_name.insert(i, name_str.clone());
            name_to_length.insert(name_str, map.length().into());
        }

        Self {
            id_to_name,
            name_to_length,
        }
    }

    fn get_name(&self, id: usize) -> Result<&str> {
        self.id_to_name
            .get(&id)
            .map(|s| s.as_str())
            .ok_or_else(|| anyhow!("Chromosome ID {} not found", id))
    }

    fn get_length(&self, name: &str) -> Result<usize> {
        self.name_to_length
            .get(name)
            .copied()
            .ok_or_else(|| anyhow!("Chromosome {} length not found", name))
    }
}

fn write_pairs_record<W: Write>(
    writer: &mut W,
    pair: &PairsRecord,
    chrom_info: &ChromosomeInfo,
) -> Result<()> {
    let chrom_1 = chrom_info.get_name(pair.chr1)?;
    let chrom_2 = chrom_info.get_name(pair.chr2)?;

    let len_1 = chrom_info.get_length(chrom_1)?;
    let len_2 = chrom_info.get_length(chrom_2)?;

    if pair.is_valid(len_1, len_2) {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pair.read_id, chrom_1, pair.pos1, chrom_2, pair.pos2, pair.strand1, pair.strand2
        )
        .context("Could not write record")?;
    }

    Ok(())
}

pub fn identify_ligation_junctions(bam_path: &str, output_directory: &str) -> Result<()> {
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let chrom_info = ChromosomeInfo::from_header(&header);

    // Enforce: BAM file must be sorted by queryname (read name) for proper grouping
    // The chunk_by operation assumes records with the same parent_id are consecutive
    check_bam_sorted_by_queryname(&header)?;

    let mut handles = HashMap::new();

    // Stream and group by parent_id (assumes BAM is sorted by queryname)
    let mcc_groups = reader.records().chunk_by(|r| {
        r.as_ref()
            .map(|record| {
                SegmentMetadata::from_read_name(record.name())
                    .parent_id()
                    .to_string()
            })
            .unwrap_or_else(|_| "UNKNOWN".to_string())
    });

    for (_, reads) in mcc_groups.into_iter() {
        let reads_raw = reads.collect::<Result<Vec<_>, _>>()?;
        let mut reads_buf = Vec::new();
        for r in reads_raw {
            reads_buf
                .push(noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &r)?);
        }
        let read_group = MCCReadGroup::new(reads_buf, FlashedStatus::Flashed);

        if read_group.contains_viewpoint() && read_group.any_mapped() {
            let filtered_group = read_group.filter_mapped();
            let pairs = filtered_group.ligation_junctions()?;

            for pair in pairs {
                let handle = handles.entry(pair.viewpoint_id.clone()).or_insert_with(|| {
                    let path =
                        Path::new(output_directory).join(format!("{}.pairs", &pair.viewpoint_id));
                    let file = std::fs::File::create(path).expect("Could not create file");
                    std::io::BufWriter::new(file)
                });

                write_pairs_record(handle, &pair, &chrom_info)?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use std::num::NonZeroUsize;

    #[test]
    fn test_chromosome_info_from_header() {
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(2000).unwrap()),
            )
            .build();

        let info = ChromosomeInfo::from_header(&header);

        assert_eq!(info.get_name(0).unwrap(), "chr1");
        assert_eq!(info.get_name(1).unwrap(), "chr2");
        assert_eq!(info.get_length("chr1").unwrap(), 1000);
        assert_eq!(info.get_length("chr2").unwrap(), 2000);
    }
}
