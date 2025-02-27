#! /usr/bin/env python

import sys
import pandas as pd
import click


def define_scan_region(
    bed6: str,
    upstream: tuple[int, int] = (0, 0),
    downstream: tuple[int, int] = (1, 500),
    itself: bool = False,
) -> pd.DataFrame:
    """
    Define the motif scan region around exons

    Args:
        bed6 (str): exon regions in BED6 format
        upstream (tuple[int, int], optional):
            upstream extension region to the exon 5' start (adjacent to 3'SS).
            Defaults to (0, 0).
        downstream (tuple[int, int], optional):
            downstream extension region to the exon 3' end (adjacent to 5'SS).
            Defaults to (1, 500).
        itself (bool, optional):
            whether include the exon sequence itself. Default to False.
    """
    exon = (
        pd.read_table(bed6, header=None, usecols=range(6))
        .set_axis(["chr", "start", "end", "name", "score", "strand"], axis=1)
        .astype({"chr": str, "start": int, "end": int, "strand": str})
        .loc[lambda x: x["strand"].eq("+") | x["strand"].eq("-")]
        .drop_duplicates()
    )

    upstream = sorted(upstream)
    if upstream[1] - upstream[0] > 0:
        upstream_scan_region = pd.concat(
            [
                exon.loc[
                    exon["strand"].eq("+")
                ].assign(  # the order end>start is important
                    end=lambda x: x["start"] - upstream[0] + 1,
                    start=lambda x: x["start"] - upstream[1],
                    name=lambda x: x["name"] + "_up",
                ),
                exon.loc[
                    exon["strand"].eq("-")
                ].assign(  # the order start>end is important
                    start=lambda x: x["end"] + upstream[0] - 1,
                    end=lambda x: x["end"] + upstream[1],
                    name=lambda x: x["name"] + "_up",
                ),
            ],
            ignore_index=True,
        )
    else:
        upstream_scan_region = None

    downstream = sorted(downstream)
    if downstream[1] - downstream[0] > 0:
        downstream_scan_region = pd.concat(
            [
                exon.loc[
                    exon["strand"].eq("+")
                ].assign(  # the order start>end is important
                    start=lambda x: x["end"] + downstream[0] - 1,
                    end=lambda x: x["end"] + downstream[1],
                    name=lambda x: x["name"] + "_dw",
                ),
                exon.loc[
                    exon["strand"].eq("-")
                ].assign(  # the order end>start is important
                    end=lambda x: x["start"] - downstream[0] + 1,
                    start=lambda x: x["start"] - downstream[1],
                    name=lambda x: x["name"] + "_dw",
                ),
            ],
            ignore_index=True,
        )
    else:
        downstream_scan_region = None

    assert not (
        upstream_scan_region is None
        and downstream_scan_region is None
        and itself is False
    )

    return (
        pd.concat(
            [upstream_scan_region, downstream_scan_region, exon if itself else None],
            ignore_index=True,
        )
        .sort_values(["chr", "start", "end", "strand"])
        .reset_index(drop=True)
    )


@click.command()
@click.option(
    "--upstream",
    "-u",
    type=(int, int),
    default=(0, 0),
    help="upstream extension region to the exon 5' start (adjacent to 3'SS). Default (0, 0)",
)
@click.option(
    "--downstream",
    "-d",
    type=(int, int),
    default=(1, 500),
    help="downstream extension region to the exon 3' end (adjacent to 5'SS). Default (1, 500)",
)
@click.option(
    "--itself",
    "-i",
    type=bool,
    default=False,
    help="whether include the exon sequence itself. Default to False",
)
@click.argument("bed6")
def cli(bed6, upstream, downstream, itself):
    """
    output the motif scan region (BED6) around exons (BED6)
    """
    define_scan_region(bed6, upstream, downstream, itself).to_csv(
        sys.stdout, sep="\t", header=False, index=False
    )


if __name__ == "__main__":
    try:
        bed6 = snakemake.input[0]
        upstream = snakemake.params["upstream"]
        downstream = snakemake.params["downstream"]
        itself = snakemake.params["itself"]
        out = snakemake.output[0]
        with open(out, "w") as O:
            define_scan_region(bed6, upstream, downstream, itself).to_csv(
                O, sep="\t", header=False, index=False
            )
    except NameError:
        cli()
