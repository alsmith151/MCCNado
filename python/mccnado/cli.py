import typer
import pathlib
from typing import Optional, List, Literal, Union
from typing_extensions import Annotated
from enum import Enum

import mccnado

app = typer.Typer()


@app.command()
def add_viewpoint_tag(bam: pathlib.Path):
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
    mccnado.annotate_bam(str(bam))

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