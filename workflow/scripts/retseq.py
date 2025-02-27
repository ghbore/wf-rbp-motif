#! /usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO
import gzip
import click


def retseq(bed4: str, reference: str) -> None:
    """
    Retrieve sequences from reference fasta based on BED regions

    Args:
        bed4 (File): BED (4 columns: chr, start, end, strand) region
        reference (File): Reference fasta file
    """
    reg = (
        pd.read_table(bed4, header=None, names=["chr", "start", "end", "strand"])
        .astype(
            {
                "chr": str,  # Important. Some chr looks like int
                "start": int,
                "end": int,
                "strand": str,
            }
        )
        .loc[lambda x: x["strand"].eq("+") | x["strand"].eq("-")]
        .drop_duplicates()
    )

    with gzip.open(reference, "rt") as IN:
        for record in SeqIO.parse(IN, "fasta"):
            reg_in_chr = reg.loc[reg["chr"].eq(record.id)].reset_index(drop=True)
            cache = [None] * reg_in_chr.shape[0]
            for i in reg_in_chr.index:
                start, end, strand = reg_in_chr.loc[i, ["start", "end", "strand"]]
                if start < 0:
                    start = 0
                if end > len(record):
                    end = len(record)
                seq = record[start:end]
                if reg_in_chr.loc[i, "strand"] == "-":
                    seq = seq.reverse_complement()
                name_galaxy_style = "_".join(
                    ["ref", record.id, str(start + 1), str(end), strand]
                )
                seq.id = name_galaxy_style
                seq.name = name_galaxy_style
                seq.description = ""
                cache[i] = seq
            SeqIO.write(cache, sys.stdout, "fasta")
            # jump out early if no more region left
            reg = reg.loc[reg["chr"].ne(record.id)]
            if reg.shape[0] == 0:
                break


@click.command()
@click.argument("bed4")
@click.argument("reference")
def cli(bed4, reference):
    retseq(bed4, reference)


if __name__ == "__main__":
    try:
        bed4 = snakemake.input[0]
        reference = snakemake.input[1]
        out_fasta = snakemake.output[0]
        with open(out_fasta, "w") as sys.stdout:
            retseq(bed4, reference)
    except NameError:
        cli()
