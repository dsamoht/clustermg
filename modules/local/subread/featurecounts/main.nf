process SUBREAD_FEATURECOUNTS {

    label 'process_low'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0' :
        'biocontainers/subread:2.1.1--h577a1d6_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(saf)

    output:
    tuple val(meta), path("*.featureCounts.tsv"), path("*.featureCounts.tsv.summary"), emit: counts
    path "versions.yml"                                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strand = [ unstranded: 0, forward: 1, reverse: 2 ][meta.strandedness]
    // Paired-end libraries are counted by fragment, once for both mates
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    """
    featureCounts \\
        -F SAF \\
        -a ${saf} \\
        -o ${prefix}.featureCounts.tsv \\
        -T ${task.cpus} \\
        -s ${strand} \\
        ${paired} \\
        ${args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(echo \$(featureCounts -v 2>&1) | sed 's/featureCounts v//')
    END_VERSIONS
    """
}
