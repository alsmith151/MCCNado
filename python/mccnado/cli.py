import typer
import pathlib
from typing import List

import mccnado

app = typer.Typer()


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
def split_viewpoint_reads():
    pass


def main():
    """
    Main function to run the CLI.
    """
    app()


if __name__ == "__main__":
    main()
