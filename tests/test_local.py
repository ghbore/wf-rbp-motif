from . import shell


@shell
def test_snakemake():
    """
    snakemake -c 2 -s workflow/Snakefile -d tests \
        results/xstreme/example/ \
        --config dna2rna=True
    """


@shell
def test_snakemake_update_report_path():
    """
    snakemake -c 2 -s workflow/Snakefile -d tests \
        reports/example/ \
        --config dna2rna=True report=reports/{prefix}
    """


@shell
def test_wdl():
    """
    cd tests/ && \
    miniwdl run -d wdl ../workflow/wf-rbp-motif.wdl \
        bed6=example.bed reference=resources/genome.fasta.gz \
        dna2rna=true \
        wf_docker=wf-rbp-motif:dev
    """
