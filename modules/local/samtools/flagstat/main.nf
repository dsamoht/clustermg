process SAMTOOLS_FLAGSTAT {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat.tsv"), emit: flagstat
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools flagstat -O tsv ${bam} > ${prefix}.flagstat.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
