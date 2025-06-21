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

- Python 3.7+
- Rust (for building from source)
- samtools (for BAM file processing)

## Usage

### Python API

```python
import mccnado

# Deduplicate FASTQ files
stats = mccnado.deduplicate_fastq(
    fastq1="input_R1.fastq.gz",
    output1="output_R1.fastq.gz",
    fastq2="input_R2.fastq.gz",  # Optional for paired-end
    output2="output_R2.fastq.gz"  # Optional for paired-end
)

print(f"Total reads: {stats['total_reads']}")
print(f"Unique reads: {stats['unique_reads']}")
print(f"Duplicate reads: {stats['duplicate_reads']}")

# Split viewpoint reads
mccnado.split_viewpoint_reads(
    bam="aligned_reads.bam",
    output="split_reads.fastq.gz"
)

# Annotate BAM file with MCC metadata
mccnado.annotate_bam(
    bam="input.bam",
    output_directory="annotated_output/"
)

# Extract ligation statistics
mccnado.extract_ligation_stats(
    bam="annotated.bam",
    stats="ligation_stats.json"
)
```

### Command Line Interface

The package also provides a command-line interface through the Python module:

```bash
# Deduplicate FASTQ files
python -m mccnado.cli deduplicate input_R1.fastq.gz output_R1.fastq.gz

# Split viewpoint reads
python -m mccnado.cli split-reads aligned_reads.bam split_reads.fastq.gz

# Annotate BAM files
python -m mccnado.cli annotate input.bam output_directory/

# Extract ligation statistics
python -m mccnado.cli ligation-stats annotated.bam stats.json
```

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