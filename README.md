# MCCNado

A high-performance Rust library with Python bindings for processing Micro-Capture-C (MCC) sequencing data.

## Overview

MCCNado is a bioinformatics tool designed for analyzing chromatin conformation capture sequencing data. It provides efficient implementations for common preprocessing tasks including FASTQ deduplication, viewpoint read splitting, BAM annotation, and ligation junction analysis.

## Features

- **FASTQ Deduplication**: Remove duplicate reads from single-end and paired-end FASTQ files
- **Viewpoint Read Splitting**: Split reads containing viewpoint sequences into constituent segments
- **BAM Annotation**: Add metadata tags to BAM files for downstream analysis
- **Ligation Junction Identification**: Extract and analyze chromatin interaction data
- **Ligation Statistics**: Generate comprehensive statistics on cis/trans interactions
- **High Performance**: Implemented in Rust with optional async processing for large datasets

## Installation

### From PyPI (recommended)
```bash
pip install mccnado
```

### From Source
```bash
git clone https://github.com/yourusername/MCCNado.git
cd MCCNado
pip install .
```

### Development Installation
```bash
git clone https://github.com/yourusername/MCCNado.git
cd MCCNado
pip install -e .
```

## Requirements

- Python 3.10+
- Rust (for building from source)
- samtools (for BAM file processing)

## Quick Start

After installation, you can immediately use the `mccnado` command:

```bash
# Deduplicate a BAM file
mccnado deduplicate-bam input.bam output.bam

# View all available commands
mccnado --help

# Get help for a specific command
mccnado deduplicate-bam --help
```

## Usage

### Tool Overview

MCCNado provides several specialized tools for MCC data processing:

#### 1. FASTQ Deduplication
Removes duplicate reads from FASTQ files by comparing sequence content and quality scores. Useful for removing PCR duplicates before alignment.

#### 2. BAM Deduplication
Removes duplicate alignments from BAM files based on genomic coordinates and alignment information. Identifies and filters PCR duplicates that have the same mapping location.

#### 3. Viewpoint Read Splitting
Splits composite reads containing viewpoint sequences into separate segments for independent analysis. Useful when reads contain both viewpoint and flanking sequence information.

#### 4. BAM Annotation
Adds MCC-specific metadata tags to BAM files, including viewpoint information, oligo coordinates, and reporter tags for classification.

#### 5. Ligation Statistics
Analyzes chromatin ligation events and generates statistics on cis/trans interactions, helping characterize the quality and type of chromatin interactions in your data.

#### 6. Ligation Junction Identification
Identifies and extracts ligation junction sequences from BAM files, useful for validating chromatin interactions and analyzing junction characteristics.

### Python API

```python
import mccnado

# 1. Deduplicate FASTQ files
# Removes duplicate reads by comparing sequences
stats = mccnado.deduplicate_fastq(
    fastq1="input_R1.fastq.gz",
    output1="output_R1.fastq.gz",
    fastq2="input_R2.fastq.gz",      # Optional for paired-end
    output2="output_R2.fastq.gz"     # Optional for paired-end
)
print(f"Total reads: {stats.total_reads}")
print(f"Unique reads: {stats.unique_reads}")
print(f"Duplicate reads: {stats.duplicate_reads}")

# 2. Deduplicate BAM files
# Removes PCR duplicates based on genomic coordinates
bam_stats = mccnado.deduplicate_bam(
    bam="aligned_reads.bam",
    output="deduplicated.bam"
)
print(f"Unique molecules: {bam_stats.unique_molecules}")
print(f"Duplicate molecules: {bam_stats.duplicate_molecules}")

# 3. Split viewpoint reads
# Separates composite reads into individual segments
mccnado.split_viewpoint_reads(
    bam="aligned_reads.bam",
    output="split_reads.bam"
)

# 4. Annotate BAM file with MCC metadata
# Adds VP (viewpoint), OC (oligo coordinates), and RT (reporter tag) tags
mccnado.annotate_bam(
    bam="input.bam",
    output="annotated.bam"
)

# 5. Extract ligation statistics
# Generates JSON report of cis/trans interactions and other statistics
mccnado.extract_ligation_stats(
    bam="annotated.bam",
    stats="ligation_stats.json"
)

# 6. Identify ligation junctions
# Extracts junction sequences and writes to output directory
mccnado.identify_ligation_junctions(
    bam="annotated.bam",
    output_directory="junctions/"
)
```

### Command Line Interface

MCCNado provides a clean, intuitive command-line interface accessible directly via the `mccnado` command after installation. The CLI uses command-line argument validation and provides helpful error messages.

#### Available Commands

```bash
# View all available commands and options
mccnado --help

# Deduplicate FASTQ files (single-end)
mccnado deduplicate-fastq input.fastq.gz output.fastq.gz

# Deduplicate FASTQ files (paired-end)
mccnado deduplicate-fastq input_R1.fastq.gz output_R1.fastq.gz \
  --fastq2 input_R2.fastq.gz --output2 output_R2.fastq.gz

# Remove PCR duplicates from BAM files
mccnado deduplicate-bam aligned_reads.bam deduplicated.bam

# Split reads containing viewpoint sequences
mccnado split-viewpoint-reads aligned_reads.bam split_reads.bam

# Annotate BAM files with MCC-specific metadata
mccnado annotate-bam input.bam annotated.bam

# Extract ligation statistics
mccnado extract-ligation-stats annotated.bam stats.json

# Identify ligation junctions
mccnado identify-ligation-junctions annotated.bam junctions/

# Get detailed help for any command
mccnado deduplicate-bam --help
mccnado deduplicate-fastq --help
```

#### CLI Features

- **Input Validation**: Automatically checks for file existence and correct file formats
- **Clear Error Messages**: Informative error reporting when issues are encountered
- **Summary Output**: Commands that deduplicate data display summary statistics
- **Help System**: Use `--help` with any command for detailed usage information

**Command Name Aliases**: Commands support both hyphenated and underscored formats (e.g., `deduplicate-bam` or `deduplicate_bam`)

## File Formats

### Input Files
- **FASTQ**: Raw sequencing reads (single-end or paired-end, gzipped or uncompressed)
- **BAM**: Aligned reads with proper headers and indexing

### Output Files
- **FASTQ**: Deduplicated reads
- **BAM**: Annotated alignment files with MCC-specific tags
- **JSON**: Ligation statistics and metadata

### BAM Tags Added by MCCNado
- `VP`: Viewpoint name
- `OC`: Oligo coordinates
- `RT`: Reporter tag (0 for capture reads, 1 for reporter reads)

## Performance

MCCNado is optimized for large-scale data processing:

- **Memory Efficient**: Streaming processing for large files
- **Parallel Processing**: Multi-threaded operations where applicable
- **Fast Hashing**: Uses xxHash for rapid duplicate detection
- **Batch Processing**: Configurable batch sizes for optimal performance

## Architecture

The package consists of several core modules:

- [`fastq_deduplicate`](src/fastq_deduplicate.rs): FASTQ deduplication logic
- [`viewpoint_read_splitter`](src/viewpoint_read_splitter.rs): Read segmentation functionality
- [`mcc_data_handler`](src/mcc_data_handler.rs): BAM annotation and processing
- [`ligation_stats`](src/ligation_stats.rs): Statistical analysis of ligation events
- [`utils`](src/utils.rs): Common utilities and data structures

## Development

### Building from Source

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Clone and build
git clone https://github.com/yourusername/MCCNado.git
cd MCCNado
cargo build --release

# Install Python package
pip install -e .
```

### Running Tests

```bash
# Rust tests
cargo test

# Python tests
python -m pytest tests/
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use MCCNado in your research, please cite:

```
[Your Citation Here]
```

## Support

For questions, issues, or feature requests, please:

1. Check the [documentation](https://github.com/yourusername/MCCNado/wiki)
2. Search existing [issues](https://github.com/yourusername/MCCNado/issues)
3. Open a new issue if needed

## Acknowledgments

- Built with [PyO3](https://pyo3.rs/) for Python-Rust interoperability
- Uses [noodles](https://github.com/zaeleus/noodles) for bioinformatics file format handling
- Powered by [tokio](https://tokio.rs/) for async operations