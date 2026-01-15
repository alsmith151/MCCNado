import pytest
import shutil
import pysam
import pathlib
import os


@pytest.fixture
def test_data_dir(tmp_path):
    return tmp_path


@pytest.fixture
def bam_file_with_duplicates(test_data_dir):
    """
    Creates a BAM file with known duplicates.
    Structure:
    - Group 1 (unique): 1 read
    - Group 2 (duplicates): 2 reads with identical coords
    - Group 3 (mixed): 2 reads with different coords
    """
    bam_path = test_data_dir / "test.bam"

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # 1. Unique read
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "A" * 10
        a.flag = (
            83  # Read 1, mapped, proper pair, reverse strand? No, just made up flags
        )
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 20
        a.cigar = ((0, 10),)
        outf.write(a)

        # 2. Duplicate reads (same group parent)
        # We need to simulate the read name structure expected by SegmentMetadata if relevant
        # But if we use simple names, we need to ensure the group parsing logic works.
        # mcc_data_handler uses: SegmentMetadata::from_read_name(record.name()).parent_id()
        # So we need proper MCC read names.

        # Example: "read_name:umi:cell_barcode:1"

        # Group 2: Duplicates
        parent_id_2 = "group2:UMI:CB"

        # Read 2a
        dup1 = pysam.AlignedSegment()
        dup1.query_name = f"{parent_id_2}:1"
        dup1.query_sequence = "C" * 10
        dup1.reference_id = 0
        dup1.reference_start = 200
        dup1.cigar = ((0, 10),)
        outf.write(dup1)

        # Read 2b (identical to 2a in coords)
        dup2 = pysam.AlignedSegment()
        dup2.query_name = f"{parent_id_2}:2"
        dup2.query_sequence = "G" * 10  # Seq diff doesn't matter for coords
        dup2.reference_id = 0
        dup2.reference_start = 200
        dup2.cigar = ((0, 10),)
        outf.write(dup2)

        # Group 3: Distinct reads in same group
        parent_id_3 = "group3:UMI:CB"

        distinct1 = pysam.AlignedSegment()
        distinct1.query_name = f"{parent_id_3}:1"
        distinct1.reference_id = 0
        distinct1.reference_start = 300
        distinct1.cigar = ((0, 10),)
        outf.write(distinct1)

        distinct2 = pysam.AlignedSegment()
        distinct2.query_name = f"{parent_id_3}:2"
        distinct2.reference_id = 0
        distinct2.reference_start = 400
        distinct2.cigar = ((0, 10),)
        outf.write(distinct2)

    pysam.index(str(bam_path))
    return bam_path
