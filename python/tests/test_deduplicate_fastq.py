import mccnado
import gzip
import pytest
from pathlib import Path


@pytest.fixture
def fastq_with_duplicates(tmp_path):
    """Create a simple FASTQ file with known duplicates."""
    fastq_file = tmp_path / "test.fastq"
    
    # Write test FASTQ with some duplicate sequences
    with open(fastq_file, "w") as f:
        # Read 1: unique
        f.write("@read1\n")
        f.write("ACGTACGTACGT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")
        
        # Read 2: duplicate of read 1 (same sequence)
        f.write("@read2\n")
        f.write("ACGTACGTACGT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")
        
        # Read 3: unique
        f.write("@read3\n")
        f.write("TGCATGCATGCA\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")
        
        # Read 4: duplicate of read 3
        f.write("@read4\n")
        f.write("TGCATGCATGCA\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")
        
        # Read 5: duplicate of read 3
        f.write("@read5\n")
        f.write("TGCATGCATGCA\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")
    
    return fastq_file


def test_deduplicate_fastq_single_end(fastq_with_duplicates, tmp_path):
    """Test single-end FASTQ deduplication."""
    output_file = tmp_path / "deduplicated.fastq"
    
    # Run deduplication
    stats = mccnado.deduplicate_fastq(
        str(fastq_with_duplicates),
        str(output_file)
    )
    
    # Verify stats
    assert stats['total_reads'] == 5
    assert stats['unique_reads'] == 2  # 1 ACGT sequence, 1 TGCA sequence
    assert stats['duplicate_reads'] == 3
    
    # Verify output file exists and has correct number of reads
    assert output_file.exists()
    
    # Count reads in output
    with open(output_file, "r") as f:
        lines = f.readlines()
        # Each read is 4 lines in FASTQ format
        num_reads = len([l for l in lines if l.startswith("@")]) 
        assert num_reads == 2


def test_deduplicate_fastq_paired_end(tmp_path):
    """Test paired-end FASTQ deduplication."""
    fastq1 = tmp_path / "test_R1.fastq"
    fastq2 = tmp_path / "test_R2.fastq"
    output1 = tmp_path / "deduplicated_R1.fastq"
    output2 = tmp_path / "deduplicated_R2.fastq"
    
    # Create paired-end FASTQ files
    with open(fastq1, "w") as f:
        f.write("@read1/1\nACGTACGT\n+\nIIIIIIII\n")
        f.write("@read2/1\nACGTACGT\n+\nIIIIIIII\n")
        f.write("@read3/1\nTGCATGCA\n+\nIIIIIIII\n")
    
    with open(fastq2, "w") as f:
        f.write("@read1/2\nTTAAGGCC\n+\nIIIIIIII\n")
        f.write("@read2/2\nTTAAGGCC\n+\nIIIIIIII\n")
        f.write("@read3/2\nCCGGAACC\n+\nIIIIIIII\n")
    
    # Run deduplication
    stats = mccnado.deduplicate_fastq(
        str(fastq1),
        str(output1),
        str(fastq2),
        str(output2)
    )
    
    # Verify stats
    assert stats['total_reads'] == 3  # Paired reads are counted as pairs
    assert stats['unique_reads'] == 2
    assert stats['duplicate_reads'] == 1
    
    # Verify both output files exist
    assert output1.exists()
    assert output2.exists()
    
    # Verify both files have same number of reads
    with open(output1, "r") as f:
        reads_r1 = len([l for l in f.readlines() if l.startswith("@")])
    with open(output2, "r") as f:
        reads_r2 = len([l for l in f.readlines() if l.startswith("@")])
    
    assert reads_r1 == reads_r2 == 2
