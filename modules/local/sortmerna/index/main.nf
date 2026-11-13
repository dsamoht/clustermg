process SORTMERNA_INDEX {

    label 'process_medium'

    tag "rRNA databases"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.7--hdbdd923_1' :
        'biocontainers/sortmerna:4.3.7--hdbdd923_1' }"

    input:
    path(fastas, stageAs: 'rrna/*')

    output:
    path "idx"         , emit: index
    path "versions.yml", emit: versions

    script:
    // Built once, then read by every SORTMERNA task, which stages the databases under the same
    // names: the index is tied to the references it was built from. A single database comes as
    // one path rather than a list.
    def refs = (fastas instanceof List ? fastas : [ fastas ]).collect { fasta -> "--ref ${fasta}" }.join(' ')
    """
    sortmerna \\
        ${refs} \\
        --workdir . \\
        --index 1 \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS
    """
}
