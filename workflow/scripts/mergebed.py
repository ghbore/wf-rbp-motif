#! /usr/bin/env python

import sys
import pandas as pd
import click


def merge1(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df["end"] > df["start"]]
    pos = (
        pd.concat(
            [
                pd.DataFrame({"pos": df["start"], "flag": 1}),
                pd.DataFrame({"pos": df["end"], "flag": -1}),
            ]
        )
        .groupby(["pos"])
        .sum()
        .sort_index()
        .cumsum()
        .rolling(2, min_periods=1)
        .apply(
            lambda x: 1
            if len(x) < 2
            else 1
            if x[0] <= 0 and x[1] > 0
            else 0
            if x[0] > 0 and x[1] <= 0
            else -1,
            raw=True,
        )
        .loc[lambda x: x["flag"] >= 0]
        .reset_index()["pos"]
    )
    return pd.DataFrame(zip(pos[0::2], pos[1::2]), columns=["start", "end"])


def mergebed(bed6: str) -> None:
    pd.read_table(bed6 if bed6 is not None else sys.stdin, header=None).iloc[
        :, :6
    ].set_axis(
        ["chr", "start", "end", "name", "score", "strand"], axis="columns"
    ).groupby(
        ["chr", "strand"]
    ).apply(
        merge1
    ).reset_index().sort_values(
        ["chr", "start", "end", "strand"]
    )[
        ["chr", "start", "end", "strand"]
    ].to_csv(
        sys.stdout, sep="\t", header=False, index=False
    )


@click.command()
@click.argument("bed6", type=click.File("r"), default=sys.stdin)
def cli(bed6):
    mergebed(bed6)


if __name__ == "__main__":
    try:
        bed6 = snakemake.input[0]
        out = snakemake.output[0]
        with open(out, "w") as sys.stdout:
            mergebed(bed6)
    except NameError:
        cli()
