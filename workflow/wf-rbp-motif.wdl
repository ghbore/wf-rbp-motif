version development


task merge_bed {
    meta {
        description: "merge user-provided BED file"
    }
    input {
        File bed6
        String docker = "ghcr.io/ghbore/wf-rbp-motif:latest"
        Int cpu = 1
        Int memory = 2
    }
    command <<<
        mergebed ~{bed6} > merged.bed4
    >>>
    output {
        File bed4 = "merged.bed4"
    }
    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GB"
    }
}


task download {
    meta {
        description: "download reference genome fasta"
    }
    input {
        String url
        Int cpu = 1
        Int memory = 2
    }
    String out = basename(url)
    command <<<
        curl --progress-bar -o ~{out} ~{url}
    >>>
    output {
        File filename = "~{out}"
    }
    runtime {
        memory: "~{memory} GB"
        cpu: cpu
    }
}


task retrieve_seq {
    meta {
        description: "retrieve sequences for motif scan"
    }
    input {
        File bed4
        File reference
        Int upstream = 1
        Int downstream = 500
        String docker = "ghcr.io/ghbore/wf-rbp-motif:latest"
        Int cpu = 1
        Int memory = 10
    }
    command <<<
        retseq --upstream ~{upstream} \
            --downstream ~{downstream} \
            ~{bed4} ~{reference} \
            > "sequence.fasta"
    >>>
    output {
        File sequence = "sequence.fasta"
    }
    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GB"
    }
}


task run_xstreme {
    meta {
        description: "Perform comprehensive motif analysis (including motif discovery) on sequences"
    }
    input {
        File sequence
        File? custom_db
        Boolean dna2rna = false
        String species_scientific_name = "Homo_sapiens"
        String docker = "memesuite/memesuite:5.5.5"
        Int cpu = 1
        Int memory = 6
    }
    String bundle_db_prefix = "/opt/meme/share/meme-*/db/motif_databases/RNA/Ray2013_rbp_"
    String bundle_db_suffix = if dna2rna then ".meme" else ".dna_encoded.meme"
    String db = select_first([
        custom_db,
        bundle_db_prefix + species_scientific_name + bundle_db_suffix
    ])
    command <<<
        xstreme \
            --p ~{sequence} \
            --oc report \
            ~{if dna2rna then "--dna2rna" else ""} \
            --seed 1 \
            --m ~{db}
    >>>
    output {
        Directory report = "report"
    }
    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GB"
    }
}


workflow RBP_Motif_Analysis {
    input {
        File bed6
        String species_scientific_name = "Homo_sapiens"
        File? reference
        String? reference_url
        Int upstream = 1
        Int downstream = 500
        File? custom_motif_db
        Boolean dna2rna = false
        String wf_docker = "ghcr.io/ghbore/wf-rbp-motif:latest"
        String mergebed_docker = "ghcr.io/ghbore/wf-rbp-motif:latest"
        String meme_docker = "ghcr.io/ghbore/wf-rbp-motif:latest"
    }
    call merge_bed {
        input:
            bed6 = bed6,
            docker = mergebed_docker
    }
    if (! defined(reference) && defined(reference_url)){
        call download as download_genome {
            input: 
                url = select_first([reference_url])
        }
    }
    call retrieve_seq {
        input:
            bed4 = merge_bed.bed4,
            reference = select_first([
                download_genome.filename,
                reference
            ]),
            upstream = upstream,
            downstream = downstream,
            docker = wf_docker
    }
    call run_xstreme {
        input:
            sequence = retrieve_seq.sequence,
            custom_db = custom_motif_db,
            dna2rna = dna2rna,
            species_scientific_name = species_scientific_name,
            docker = meme_docker
    }
    output {
        Directory report = run_xstreme.report
    }
}