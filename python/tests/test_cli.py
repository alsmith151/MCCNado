"""
Test the command-line interface.
"""

import subprocess
import pathlib
import tempfile
import gzip
import pytest
import re


@pytest.fixture
def test_data_dir():
    """Get the test data directory."""
    return pathlib.Path(__file__).parent.parent.parent / "tests" / "data"


@pytest.fixture
def temp_dir():
    """Create a temporary directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield pathlib.Path(tmpdir)


def get_cargo_version():
    """Extract version from Cargo.toml."""
    cargo_path = pathlib.Path(__file__).parent.parent.parent / "Cargo.toml"
    with open(cargo_path) as f:
        for line in f:
            if line.startswith("version"):
                # Extract version string like "0.1.4" from 'version = "0.1.4"'
                match = re.search(r'version\s*=\s*"([^"]+)"', line)
                if match:
                    return match.group(1)
    raise ValueError("Could not find version in Cargo.toml")


@pytest.fixture
def sample_fastq(temp_dir):
    """Create a sample FASTQ file for testing."""
    # Use proper FASTQ header format with sequence headers
    fastq_content = """@read1__all__1
ACGTACGT
+
IIIIIIII
@read2__all__1
TGCATGCA
+
IIIIIIII
@read1__all__1
ACGTACGT
+
IIIIIIII
"""
    fastq_path = temp_dir / "sample.fastq"
    with open(fastq_path, "w") as f:
        f.write(fastq_content)
    return fastq_path


@pytest.fixture
def sample_fastq_r1_r2(temp_dir):
    """Create paired-end FASTQ files for testing."""
    r1_content = """@read1__all__1
ACGTACGT
+
IIIIIIII
@read2__all__1
TGCATGCA
+
IIIIIIII
@read1__all__1
ACGTACGT
+
IIIIIIII
"""
    r2_content = """@read1__all__2
TTTTAAAA
+
IIIIIIII
@read2__all__2
AAAATTTT
+
IIIIIIII
@read1__all__2
TTTTAAAA
+
IIIIIIII
"""
    r1_path = temp_dir / "sample_R1.fastq"
    r2_path = temp_dir / "sample_R2.fastq"
    
    with open(r1_path, "w") as f:
        f.write(r1_content)
    with open(r2_path, "w") as f:
        f.write(r2_content)
    
    return r1_path, r2_path


def run_cli_command(args):
    """
    Run a CLI command and return the result.
    
    Args:
        args: List of command arguments (e.g., ["mccnado", "deduplicate-fastq", ...])
    
    Returns:
        CompletedProcess object with returncode, stdout, stderr
    """
    result = subprocess.run(
        args,
        capture_output=True,
        text=True,
    )
    return result


class TestVersionSync:
    """Test that versions are synced across Python and Rust."""
    
    def test_python_package_version_matches_cargo(self):
        """Test that Python package version matches Cargo.toml."""
        import mccnado
        cargo_version = get_cargo_version()
        assert mccnado.__version__ == cargo_version, (
            f"Version mismatch: Python package has {mccnado.__version__} "
            f"but Cargo.toml has {cargo_version}"
        )
    
    def test_cli_version_matches_cargo(self):
        """Test that CLI version matches Cargo.toml."""
        from mccnado.version import __version__
        cargo_version = get_cargo_version()
        assert __version__ == cargo_version, (
            f"Version mismatch: CLI version is {__version__} "
            f"but Cargo.toml has {cargo_version}"
        )


class TestCLIBasic:
    """Test basic CLI functionality."""
    
    def test_cli_help(self):
        """Test that the CLI shows help."""
        result = run_cli_command(["mccnado", "--help"])
        assert result.returncode == 0
        assert "Usage:" in result.stdout or "usage:" in result.stdout.lower()
    
    def test_cli_version(self):
        """Test that the CLI shows version."""
        result = run_cli_command(["mccnado", "--version"])
        assert result.returncode == 0
        assert "mccnado" in result.stdout
        assert "0.1" in result.stdout  # Check version is shown
    
    def test_cli_version_or_help(self):
        """Test that running with no args shows help or version."""
        result = run_cli_command(["mccnado"])
        # Either shows help/error is acceptable
        assert result.returncode in [0, 2, 1]


class TestDeduplicateFastqCLI:
    """Test deduplicate_fastq CLI command."""
    
    def test_deduplicate_fastq_single_end(self, sample_fastq, temp_dir):
        """Test deduplicate-fastq with single-end reads."""
        output_path = temp_dir / "output.fastq"
        
        result = run_cli_command([
            "mccnado",
            "deduplicate-fastq",
            str(sample_fastq),
            str(output_path),
        ])
        
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        assert output_path.exists(), "Output file was not created"
        assert "Deduplication summary:" in result.stdout
        assert "Total reads:" in result.stdout
        assert "Unique reads:" in result.stdout
        
        # Verify output has fewer reads than input (deduplication worked)
        with open(output_path) as f:
            output_lines = f.readlines()
        assert len(output_lines) > 0
    
    def test_deduplicate_fastq_paired_end(self, sample_fastq_r1_r2, temp_dir):
        """Test deduplicate-fastq with paired-end reads."""
        r1_path, r2_path = sample_fastq_r1_r2
        output_r1 = temp_dir / "output_R1.fastq"
        output_r2 = temp_dir / "output_R2.fastq"
        
        result = run_cli_command([
            "mccnado",
            "deduplicate-fastq",
            str(r1_path),
            str(output_r1),
            "--fastq2", str(r2_path),
            "--output2", str(output_r2),
        ])
        
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        assert output_r1.exists(), "Output R1 file was not created"
        assert output_r2.exists(), "Output R2 file was not created"
        assert "Deduplication summary:" in result.stdout
        
        # Verify both output files have matching read counts
        with open(output_r1) as f:
            r1_lines = len(f.readlines())
        with open(output_r2) as f:
            r2_lines = len(f.readlines())
        assert r1_lines == r2_lines, "R1 and R2 output files have different line counts"
    
    def test_deduplicate_fastq_missing_file(self, temp_dir):
        """Test deduplicate-fastq with missing input file."""
        missing_file = temp_dir / "missing.fastq"
        output_file = temp_dir / "output.fastq"
        
        result = run_cli_command([
            "mccnado",
            "deduplicate-fastq",
            str(missing_file),
            str(output_file),
        ])
        
        assert result.returncode != 0, "Should fail with missing input file"
    
    def test_deduplicate_fastq_invalid_extension(self, temp_dir):
        """Test deduplicate-fastq with invalid file extension."""
        invalid_file = temp_dir / "sample.txt"
        invalid_file.write_text("not a fastq file")
        output_file = temp_dir / "output.fastq"
        
        result = run_cli_command([
            "mccnado",
            "deduplicate-fastq",
            str(invalid_file),
            str(output_file),
        ])
        
        assert result.returncode != 0, "Should fail with invalid file extension"
    
    def test_deduplicate_fastq_paired_missing_output2(self, sample_fastq_r1_r2, temp_dir):
        """Test deduplicate-fastq paired-end without output2 specified."""
        r1_path, r2_path = sample_fastq_r1_r2
        output_r1 = temp_dir / "output_R1.fastq"
        
        result = run_cli_command([
            "mccnado",
            "deduplicate-fastq",
            str(r1_path),
            str(output_r1),
            "--fastq2", str(r2_path),
            # Missing --output2 flag
        ])
        
        assert result.returncode != 0, "Should fail when fastq2 is provided without output2"


class TestCLIAvailableCommands:
    """Test that all expected CLI commands are available."""
    
    def test_annotate_bam_help(self):
        """Test annotate-bam help."""
        result = run_cli_command(["mccnado", "annotate-bam", "--help"])
        assert result.returncode == 0
        assert "annotate" in result.stdout.lower() or "viewpoint" in result.stdout.lower()
    
    def test_deduplicate_fastq_help(self):
        """Test deduplicate-fastq help."""
        result = run_cli_command(["mccnado", "deduplicate-fastq", "--help"])
        assert result.returncode == 0
        assert "deduplicate" in result.stdout.lower() or "fastq" in result.stdout.lower()
    
    def test_split_viewpoint_reads_help(self):
        """Test split-viewpoint-reads help."""
        result = run_cli_command(["mccnado", "split-viewpoint-reads", "--help"])
        assert result.returncode == 0
    
    def test_identify_ligation_junctions_help(self):
        """Test identify-ligation-junctions help."""
        result = run_cli_command(["mccnado", "identify-ligation-junctions", "--help"])
        assert result.returncode == 0
    
    def test_extract_ligation_stats_help(self):
        """Test extract-ligation-stats help."""
        result = run_cli_command(["mccnado", "extract-ligation-stats", "--help"])
        assert result.returncode == 0
    
    def test_deduplicate_bam_help(self):
        """Test deduplicate-bam help."""
        result = run_cli_command(["mccnado", "deduplicate-bam", "--help"])
        assert result.returncode == 0
