from . import shell


@shell
def test_snakemake():
    """
    docker run --rm \
        -v $PWD:$PWD -w $PWD \
        wf-rbp-motif:dev \
            snakemake -c 2 -s /wf/workflow/Snakefile \
                -d tests results/xstreme/example-docker/ \
                --config dna2rna=True report=results/xstreme/{prefix}-docker/
    """


@shell
def test_wdl():
    """
    cd tests/ && \
    docker run --rm \
        -v $PWD:$PWD -w $PWD \
        wf-rbp-motif:dev \
            miniwdl run -d wdl /wf/workflow/wf-rbp-motif.wdl \
                bed6=example.bed reference=resources/genome.fasta.gz \
                dna2rna=true
    """
