import pandas as pd
import gzip


stat_url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05820-3/MediaObjects/41586_2023_5820_MOESM5_ESM.xlsx"
gff_url = "https://ftp.ensembl.org/pub/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh38.87.gff3.gz"
prefix = "Jbara-nature-2023-extfig2a"


configfile: "jbara.yaml"


rule all:
    input:
        f"results/xstreme/{prefix}/",


rule download_as_stat:
    output:
        f"resources/{prefix}.xlsx",
    params:
        url=stat_url,
    shell:
        """
        curl --progress-bar -o {output} {params.url}
        """


use rule download_as_stat as download_gff with:
    output:
        "resources/hg38.gff.gz",
    params:
        url=gff_url,


rule extract_target_exons:
    input:
        ancient(rules.download_as_stat.output),
    output:
        f"results/{prefix}.target_exon",
    run:
        pd.read_excel(
            input[0],
            usecols=[
                "Target Exon",
                        "Reference Transcript",
                        "Î”PSI (%)",
                        "T-test p-value",
                        "FDR (BH)",
                    ],
                    na_values="na",
                ).set_axis(
                ["exon", "transcript", "delta_psi", "pval", "padj"], axis="columns"
        ).assign(
            transcript=lambda x: x["transcript"].str.replace(r"^.+\.", "", regex=True)
        ).loc[
            # Jbara et al. applied threshold on "nominal P value",
            # but here, prefer to "adjusted P value (BH)"
            lambda x: x["delta_psi"].abs().gt(10)
            & x["padj"].lt(0.05)
        ].to_csv(
            output[0],
            columns=["exon", "transcript"],
            header=True,
            index=False,
            sep="\t",
        )


def load_gff(file: str) -> pd.DataFrame:
    """load GFF3 file as DataFrame"""
    cache = dict()
    with gzip.open(file, "rt") as GFF:
        for line in GFF:
            if "#" != line[0]:
                cols = line.split("\t")
                anno = dict(term.split("=") for term in cols[8].split(";"))
                cache[anno.get("transcript_id", "")] = {
                    "transcript": anno.get("transcript_id", ""),
                    "strand": cols[6],
                }
    if "" in cache:
        del cache[""]
    return pd.DataFrame(cache.values()).drop_duplicates()


rule add_strand_info:
    input:
        rules.extract_target_exons.output,
        gff=ancient(rules.download_gff.output),
    output:
        f"results/{prefix}.bed",
    run:
        df = (
            pd.read_table(input[0])
            .merge(
                load_gff(str(input["gff"])).assign(
                    transcript=lambda x: x["transcript"].str.replace(
                        r"\.[0-9]+$", "", regex=True
                    )
                ),
                how="left",
            )
            .loc[lambda x: x["strand"].eq("+") | x["strand"].eq("-")]
        )
        df[["chr", "start", "end"]] = df["exon"].str.split(
            r":|-", expand=True, regex=True
        )

        # dump
        df[["chr", "start", "end", "strand"]].astype({"start": int, "end": int}).assign(
            start=lambda x: x["start"] - 1, name=".", score=0
        )[["chr", "start", "end", "name", "score", "strand"]].sort_values(
            ["chr", "start", "end", "strand"]
        ).drop_duplicates(
            ["chr", "start", "end", "strand"]
        ).to_csv(
            output[0], header=False, index=False, sep="\t"
        )


module meme:
    snakefile:
        github("ghbore/wf-rbp-motif", path="workflow/Snakefile", tag="main")
    config:
        config


use rule * from meme as meme_*
