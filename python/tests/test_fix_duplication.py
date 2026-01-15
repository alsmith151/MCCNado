import mccnado
import pysam
import pytest
import os


def test_deduplicate_cleans_rt_duplicates(tmp_path):
    bam_path = tmp_path / "rt_duplicate_test.bam"
    output_bam = tmp_path / "rt_deduplicated.bam"

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    # Create BAM with a read that has both RT=0 and RT=1 entries (simulating the bug)
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Parent ID structure: parent__viewpoint__pos__readnum__flashed
        # Using a valid name structure
        read_name = "parent1__vp1__start__1__1:1"

        # Record 1: RT=0 (Capture)
        a = pysam.AlignedSegment()
        a.query_name = read_name
        a.query_sequence = "A" * 10
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.cigar = ((0, 10),)
        a.set_tag("RT", 0)
        outf.write(a)  # Write first copy

        # Record 2: RT=1 (Reporter) - Identical otherwise
        b = pysam.AlignedSegment()
        b.query_name = read_name
        b.query_sequence = "A" * 10
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 100
        b.cigar = ((0, 10),)
        b.set_tag("RT", 1)
        outf.write(b)  # Write second copy

    pysam.index(str(bam_path))

    # Run deduplication
    stats = mccnado.deduplicate_bam(str(bam_path), str(output_bam))

    # Verify stats
    # Should see 1 unique molecule (since coords are identical and deduped in key)
    assert stats.total_molecules == 1, (
        f"Expected 1 total molecule group, got {stats.total_molecules}"
    )
    assert stats.unique_molecules == 1

    # Verify output content
    with pysam.AlignmentFile(output_bam, "rb") as bam:
        reads = list(bam)
        assert len(reads) == 1, "Should have filtered down to 1 record"

        read = reads[0]
        rt_tag = read.get_tag("RT")
        assert rt_tag == 1, f"Expected RT=1 (Reporter), got RT={rt_tag}"


def test_deduplicate_handles_valid_duplicates(tmp_path):
    # Verify we still deduplicate distinct PCR duplicates
    bam_path = tmp_path / "pcr_dup_test.bam"
    output_bam = tmp_path / "pcr_deduplicated.bam"

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Molecule 1
        name1 = "mol1__vp__start__1__1:1"
        a = pysam.AlignedSegment()
        a.query_name = name1
        a.reference_id = 0
        a.reference_start = 100
        a.set_tag("RT", 1)
        outf.write(a)

        # Molecule 2 (PCR Duplicate of 1)
        name2 = "mol2__vp__start__1__1:1"  # Different parent ID
        b = pysam.AlignedSegment()
        b.query_name = name2
        b.reference_id = 0
        b.reference_start = 100  # Same start
        b.set_tag("RT", 1)
        outf.write(b)

    pysam.index(str(bam_path))

    stats = mccnado.deduplicate_bam(str(bam_path), str(output_bam))

    assert stats.unique_molecules == 1
    assert stats.duplicate_molecules == 1
