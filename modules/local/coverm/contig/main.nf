process COVERM_CONTIG {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hcb7b614_4' :
        'biocontainers/coverm:0.7.0--hcb7b614_4' }"

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("coverm_contig.tsv"), emit: stats
    path "versions.yml"                       , emit: versions

    script:
    // Every contig of the reference, in every DNA library: the depth its genes are given. A
    // column is named after its BAM, <library> Mean and <library> Trimmed Mean.
    def args = task.ext.args ?: ''
    """
    coverm contig \\
        ${args} \\
        --bam-files ${bams} \\
        --threads ${task.cpus} \\
        --output-file coverm_contig.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version | sed 's/coverm //')
    END_VERSIONS
    """
}
