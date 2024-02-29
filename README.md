# Workflow for Preforming Motif Analysis of RNA Binding Proteins

This repo contains [a WDL workflow](workflow/wf-rbp-motif.wdl) and [a Snakemake workflow](workflow/Snakefile) to perform comprehensive RNA protein binding motif analysis (including motif discovery) on sequences, by wrapping [xstreme](#1) method in [the MEME suite](https://meme-suite.org/meme/doc/xstreme.html?man_type=web).

The primary inputs include
  1. a BED6+ file containing the target regions and 
  2. a reference FASTA file

Users has the option to provide a pre-defined motif database, either in the form of [a custom MEME database](https://meme-suite.org/meme/doc/meme-format.html) or by using the MEME suite bundled [Ray2013 RBP motif](#2) database.

For the [Snakemake workflow](workflow/Snakefile), the configuration YAML looks like:

```YAML
__use_yte__: true

__variables__:
  scientific_name:
    human: Homo_sapiens
    mouse: Mus_musculus

species: human
species_scientific_name: ?scientific_name[this["species"]]
genome:
  url: "https://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
  local: "resources/hg38.fasta.gz"

scan_region:
  upstream: 1
  downstream: 500

dna2rna: False
```

  1. *species_scientific_name* specifies the scientific name of the species, which is used to generate the path to the bundled Ray2013 RBP motif database for that species.
  2. *genome* configures both the remote url and local path of the reference genome FASTA file.
  3. *scan_region* defines how much to expand the target regions for analysis.
  4. *dna2rna* is used in *XSTREME* program, indicating that the input DNA sequences will be treated as RNA (treating T as representing U), and the output motifs will use the RNA alphabet.

[The included demo](./jbara.smk) intends to replicate [the extend figure 2a](https://www.nature.com/articles/s41586-023-05820-3/figures/6) from [Jbara2023](#3).

To run this demo, 

```bash
docker run --rm \
  -v $PWD:$PWD -w $PWD \
  ghcr.io/ghbore/wf-rbp-motif:latest \
    snakemake -c all -s /opt/wf/jbara.smk \
      --configfile /opt/wf/jbara.yaml
```


For the [WDL workflow](workflow/wf-rbp-motif.wdl), the possible inputs look like:

```WDL
File bed6
String species_scientific_name = "Homo_sapiens"
File? reference
String? reference_url
Int upstream = 1
Int downstream = 500
File? custom_motif_db
Boolean dna2rna = false
String wf_docker = "ghcr.io/ghbore/wf-rbp-motif:latest"
String bedtools_docker = "biocontainers/bedtools:2.25.0"
String meme_docker = "memesuite/memesuite:5.5.5"
```

Most parameters align with those in the Snakemake workflow. Similarly, [a demo input json](./jbara-input.json) is provided to replicate [the extend figure 2a](https://www.nature.com/articles/s41586-023-05820-3/figures/6) from [Jbara2023](#3), using the intermediate files generated by [the Snakemake demo](./jbara.smk).

To run this demo, 

```bash
tree  # dependent files
# resources/
# ├── Jbara-nature-2023-extfig2a.xlsx
# ├── hg38.fasta.gz
# └── hg38.gff.gz
# results/
# ├── Jbara-nature-2023-extfig2a.bed
# ├── Jbara-nature-2023-extfig2a.fasta
# ├── Jbara-nature-2023-extfig2a.merged.bed
# ├── Jbara-nature-2023-extfig2a.target_exon
# ├── xstreme
# │   └── Jbara-nature-2023-extfig2a

docker run --rm \
  -v $PWD:$PWD -w $PWD \
  ghcr.io/ghbore/wf-rbp-motif \
    miniwdl run \
      -i /opt/wf/jbara-input.json \
      -o jbara-output.json \
      /opt/wf/workflow/wf-rbp-motif.wdl
```


## Reference

<a id="1">1</a>. Charles E. Grant and Timothy L. Bailey, "XSTREME: comprehensive motif analysis of biological sequence datasets", BioRxiv, September 3, 2021.

<a id="2">2</a>. Ray D, Kazan H, Cook KB, Weirauch MT, Najafabadi HS, Li X, Gueroussov S, Albu M, Zheng H, Yang A, Na H, Irimia M, Matzat LH, Dale RK, Smith SA, Yarosh CA, Kelly SM, Nabet B, Mecenas D, Li W, Laishram RS, Qiao M, Lipshitz HD, Piano F, Corbett AH, Carstens RP, Frey BJ, Anderson RA, Lynch KW, Penalva LO, Lei EP, Fraser AG, Blencowe BJ, Morris QD, Hughes TR. A compendium of RNA-binding motifs for decoding gene regulation. Nature. 2013 Jul 11;499(7457):172-7. doi: 10.1038/nature12311IF: 64.8 Q1 . PMID: 23846655; PMCID: PMC3929597.

<a id="3">3.</a> Jbara A, Lin KT, Stossel C, Siegfried Z, Shqerat H, Amar-Schwartz A, Elyada E, Mogilevsky M, Raitses-Gurevich M, Johnson JL, Yaron TM, Ovadia O, Jang GH, Danan-Gotthold M, Cantley LC, Levanon EY, Gallinger S, Krainer AR, Golan T, Karni R. RBFOX2 modulates a metastatic signature of alternative splicing in pancreatic cancer. Nature. 2023 May;617(7959):147-153. doi: 10.1038/s41586-023-05820-3. Epub 2023 Mar 22. PMID: 36949200; PMCID: PMC10156590.