process MMSEQS_CLUSTER {

    label 'process_high'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0' :
        'biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0' }"

    input:
    tuple val(meta), path(proteins, stageAs: 'input/*')

    output:
    tuple val(meta), path("*_cluster.tsv")  , emit: clusters
    tuple val(meta), path("*_rep_seq.fasta"), emit: rep_seq
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The member sequences are already in the inputs: `_all_seqs.fasta` is a copy of every
    // input sequence, removed rather than left in the work directory
    """
    mmseqs easy-cluster \\
        input/* \\
        ${prefix} \\
        tmp \\
        --threads ${task.cpus} \\
        ${args}

    rm -rf tmp ${prefix}_all_seqs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """
}
