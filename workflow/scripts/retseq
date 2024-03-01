#! /usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO
import gzip
import click


def retseq(bed4, reference, upstream=1, downstream=500):
    """
    Retrieve sequences from reference fasta based on BED regions

    Args:
        bed4 (File): BED (4 columns: chr, start, end, strand) region
        reference (File): Reference fasta file

    Options:
        upstream (Int): # of bases extended to the upstream
        downstream (Int): # of bases extended to the downstream
    """
    reg = (
        pd.read_table(bed4, header=None, names=["chr", "start", "end", "strand"])
        .astype({"chr": str, "start": int, "end": int, "strand": str})
        .loc[lambda x: x["strand"].eq("+") | x["strand"].eq("-")]
        .drop_duplicates()
    )

    with gzip.open(reference, "rt") as IN:
        for record in SeqIO.parse(IN, "fasta"):
            reg_in_chr = reg.loc[reg["chr"].eq(record.id)].reset_index(drop=True)
            cache = [None] * reg_in_chr.shape[0]
            for i in reg_in_chr.index:
                if reg_in_chr.loc[i, "strand"] == "+":
                    start, end = reg_in_chr.loc[i, ["start", "start"]] + [
                        upstream - 1,
                        downstream,
                    ]
                    seq = record[start:end]
                else:
                    start, end = reg_in_chr.loc[i, ["end", "end"]] - [
                        downstream,
                        upstream - 1,
                    ]
                    seq = record[start:end].reverse_complement()
                name_galaxy_style = "_".join(
                    [
                        "ref",
                        record.id,
                        str(start),
                        str(end),
                        reg_in_chr.loc[i, "strand"],
                    ]
                )
                seq.id = name_galaxy_style
                seq.name = name_galaxy_style
                seq.description = ""
                cache[i] = seq
            SeqIO.write(cache, sys.stdout, "fasta")


@click.command()
@click.option("--upstream", "-u", default=1)
@click.option("--downstream", "-d", default=500)
@click.argument("bed4")
@click.argument("reference")
def cli(bed4, reference, upstream, downstream):
    retseq(bed4, reference, upstream, downstream)


if __name__ == "__main__":
    try:
        bed4 = snakemake.input[0]
        reference = snakemake.input[1]
        upstream = snakemake.params["upstream"]
        downstream = snakemake.params["downstream"]
        out_fasta = snakemake.output[0]
        with open(out_fasta, "w") as sys.stdout:
            retseq(bed4, reference, upstream, downstream)
    except NameError:
        cli()
