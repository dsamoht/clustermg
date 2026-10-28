process MMSEQS_SEARCH {

    label 'process_high'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0' :
        'biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0' }"

    input:
    tuple val(meta), path(queries, stageAs: 'query/*')
    path(targets, stageAs: 'target/*')

    output:
    tuple val(meta), path("*.m8"), emit: hits
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mmseqs easy-search \\
        query/* \\
        target/* \\
        ${prefix}.m8 \\
        tmp \\
        --threads ${task.cpus} \\
        --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits \\
        ${args}

    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """
}
