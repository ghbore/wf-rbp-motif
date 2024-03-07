from . import shell, docker


@shell
def test_snakemake(docker):
    """
    docker run --rm \
        -v $PWD:$PWD -w $PWD \
        wf-rbp-motif:dev \
            snakemake -c 2 -s /wf/workflow/Snakefile \
                -d tests results/xstreme/example-docker/ \
                --config dna2rna=True report=results/xstreme/{prefix}-docker/
    """


@shell
def test_wdl(docker):
    """
    cd tests/ && \
    docker run --rm \
        -v $PWD:$PWD -w $PWD \
        wf-rbp-motif:dev \
            miniwdl run -d wdl /wf/workflow/wf-rbp-motif.wdl \
                bed6=example.bed reference=resources/genome.fasta.gz \
                dna2rna=true
    """
