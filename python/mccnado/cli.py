import typer
import pathlib
from typing import List, Optional

import mccnado
from mccnado.version import __version__

app = typer.Typer()


def version_callback(value: bool):
    """Callback for --version flag."""
    if value:
        print(f"mccnado {__version__}")
        raise typer.Exit()


@app.callback()
def main_callback(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit",
    ),
):
    """MCCNado: Rust-based tools for processing Micro-Capture-C data."""
    pass


def check_bam_sorted_by_queryname(bam_path: str):
    import pysam

    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        sorted_order = bam_file.header.get("HD", {}).get("SO")

        if sorted_order != "queryname":
            raise ValueError(
                f"The BAM file {bam_path} must be sorted by queryname (read name), currently {sorted_order}. Please sort the BAM file accordingly."
            )



@app.command()
def annotate_bam(bam: pathlib.Path, output: pathlib.Path):
    """
    Add a viewpoint tag to the BAM file.
    """
    # Check if the BAM file exists
    if not bam.exists():
        raise FileNotFoundError(f"The file {bam} does not exist.")

    # Check if the file is a BAM file
    if bam.suffix != ".bam":
        raise ValueError(f"The file {bam} is not a BAM file.")
    
    # Check that the bam files is queryname sorted
    check_bam_sorted_by_queryname(str(bam))

    # Add the viewpoint tag to the BAM file
    mccnado.annotate_bam(str(bam), str(output))


@app.command()
def extract_ligation_stats(bam: pathlib.Path, stats: pathlib.Path):
    """
    Extract ligation statistics from the BAM file.
    """
    # Check if the BAM file exists
    if not bam.exists():
        raise FileNotFoundError(f"The file {bam} does not exist.")

    # Check if the file is a BAM file
    if bam.suffix != ".bam":
        raise ValueError(f"The file {bam} is not a BAM file.")

    # Extract ligation statistics from the BAM file
    mccnado.extract_ligation_stats(str(bam), str(stats))


@app.command()
def identify_ligation_junctions(bam: pathlib.Path, outdir: pathlib.Path):
    """
    Identify ligation junctions from the BAM file.
    """
    # Check if the BAM file exists
    if not bam.exists():
        raise FileNotFoundError(f"The file {bam} does not exist.")
    # Check if the file is a BAM file
    if bam.suffix != ".bam":
        raise ValueError(f"The file {bam} is not a BAM file.")
    # Check if the output directory exists if not, create it
    if not outdir.exists():
        outdir.mkdir(parents=True)
    
    # Check that the bam files is queryname sorted
    check_bam_sorted_by_queryname(str(bam))

    # Identify ligation junctions from the BAM file
    mccnado.identify_ligation_junctions(str(bam), str(outdir))


@app.command()
def deduplicate_bam(bam: pathlib.Path, output: pathlib.Path):
    """
    Remove duplicate molecules from a BAM file based on segment coordinates.
    """
    # Check if the BAM file exists
    if not bam.exists():
        raise FileNotFoundError(f"The file {bam} does not exist.")

    # Check if the file is a BAM file
    if bam.suffix != ".bam":
        raise ValueError(f"The file {bam} is not a BAM file.")

    # Check that the bam files is queryname sorted
    check_bam_sorted_by_queryname(str(bam))

    # Deduplicate the BAM file
    stats = mccnado.deduplicate_bam(str(bam), str(output))
    print("Deduplication summary:")
    print(f"  Total molecules:     {stats.total_molecules}")
    print(f"  Unique molecules:    {stats.unique_molecules}")
    print(f"  Duplicate molecules: {stats.duplicate_molecules}")


@app.command()
def combine_ligation_junction_coolers(
    clrs: List[pathlib.Path],
    outfile: pathlib.Path,
):
    """
    Combine ligation junctions from multiple Cooler formatted files into a single file.
    """
    from .storage import CoolerBinsLinker, CoolerMerger

    # Check if the Cooler files exist
    for clr in clrs:
        if not clr.exists():
            raise FileNotFoundError(f"The file {clr} does not exist.")
        # Check if the file is a Cooler file
        if clr.suffix not in [".cool", ".mcool"]:
            raise ValueError(f"The file {clr} is not a Cooler file.")

    # Combine the Cooler files -- TODO: allow for names to be passed in
    clr_merger = CoolerMerger(clrs, outfile)
    clr_merger.merge()
    # Check if the output file exists
    if not outfile.exists():
        raise FileNotFoundError(
            f"The file {outfile} does not exist. Error merging files."
        )
    # Link the bins in the Cooler file to save space
    clr_bins_linker = CoolerBinsLinker(outfile)
    clr_bins_linker.link_bins()


@app.command()
def split_viewpoint_reads(bam: pathlib.Path, output: pathlib.Path):
    """
    Split reads containing viewpoint sequences into constituent segments.
    """
    # Check if the BAM file exists
    if not bam.exists():
        raise FileNotFoundError(f"The file {bam} does not exist.")
    
    # Check if the file is a BAM file
    if bam.suffix != ".bam":
        raise ValueError(f"The file {bam} is not a BAM file.")
    
    # Split viewpoint reads from the BAM file
    mccnado.split_viewpoint_reads(str(bam), str(output))


@app.command()
def deduplicate_fastq(
    fastq1: pathlib.Path,
    output1: pathlib.Path,
    fastq2: pathlib.Path = typer.Option(None, "--fastq2", help="Optional R2 FASTQ file"),
    output2: pathlib.Path = typer.Option(None, "--output2", help="Optional output R2 FASTQ file"),
):
    """
    Remove duplicate reads from FASTQ files.
    """
    # Check if the FASTQ files exist
    if not fastq1.exists():
        raise FileNotFoundError(f"The file {fastq1} does not exist.")
    
    # Validate file extensions
    valid_extensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    if fastq1.suffix not in valid_extensions and not any(str(fastq1).endswith(ext) for ext in valid_extensions):
        raise ValueError(f"The file {fastq1} does not appear to be a FASTQ file.")
    
    if fastq2 is not None:
        if not fastq2.exists():
            raise FileNotFoundError(f"The file {fastq2} does not exist.")
        if output2 is None:
            raise ValueError("output2 is required when fastq2 is provided")
        
        # Deduplicate paired-end FASTQ files
        stats = mccnado.deduplicate_fastq(str(fastq1), str(output1), str(fastq2), str(output2))
    else:
        # Deduplicate single-end FASTQ file
        stats = mccnado.deduplicate_fastq(str(fastq1), str(output1))
    
    # Print deduplication summary
    print("Deduplication summary:")
    print(f"  Total reads:     {stats['total_reads']}")
    print(f"  Unique reads:    {stats['unique_reads']}")
    print(f"  Duplicate reads: {stats['duplicate_reads']}")


def main():
    """
    Main function to run the CLI.
    """
    app()


if __name__ == "__main__":
    main()
