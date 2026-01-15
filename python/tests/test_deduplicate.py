import mccnado
import pysam
import os


def test_deduplicate_bam(bam_file_with_duplicates, tmp_path):
    output_bam = tmp_path / "deduplicated.bam"

    stats = mccnado.deduplicate_bam(str(bam_file_with_duplicates), str(output_bam))

    # Check stats
    assert stats.total_molecules == 5
    assert stats.unique_molecules == 4  # 1 (Group 1) + 1 (Group 2) + 2 (Group 3)
    assert stats.duplicate_molecules == 1  # 1 from Group 2

    # Check output BAM content
    assert output_bam.exists()

    with pysam.AlignmentFile(output_bam, "rb") as bam:
        reads = list(bam)
        assert len(reads) == 4

        # Verify read names present
        read_names = sorted([r.query_name for r in reads])
        expected_names = sorted(
            [
                "read1",
                "group2:UMI:CB:1",  # Only first one kept
                "group3:UMI:CB:1",
                "group3:UMI:CB:2",
            ]
        )
        assert read_names == expected_names
