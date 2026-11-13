process COVERM_GENOME {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hcb7b614_4' :
        'biocontainers/coverm:0.7.0--hcb7b614_4' }"

    input:
    tuple val(meta), path(bams), path(bais)
    path(mags, stageAs: 'mags/*')

    output:
    tuple val(meta), path("coverm_genome.tsv"), emit: stats
    path "versions.yml"                       , emit: versions

    script:
    // Every MAG, on its own contigs only: the profile of each version of a genome. The reads
    // of a whole dRep cluster are added up from the contig counts instead.
    def args = task.ext.args ?: ''
    """
    coverm genome \\
        ${args} \\
        --genome-fasta-directory mags \\
        --genome-fasta-extension fa \\
        --bam-files ${bams} \\
        --threads ${task.cpus} \\
        --output-file coverm_genome.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version | sed 's/coverm //')
    END_VERSIONS
    """
}
