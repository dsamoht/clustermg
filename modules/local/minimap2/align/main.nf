process MINIMAP2_ALIGN {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:eb7d926a9d7da1ab8e5427cd947dbbb8795eabab-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:eb7d926a9d7da1ab8e5427cd947dbbb8795eabab-0' }"

    input:
    tuple val(meta), path(reads), path(reference)
    val bam_format

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam, optional: true
    tuple val(meta), path("*.paf")                   , emit: paf, optional: true
    path "versions.yml"                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A read library brings its preset in its meta; other alignments give theirs in ext.args
    def preset = meta.preset ? "-x ${meta.preset}" : ''
    if (bam_format) {
        // samtools sort takes up to 4 of the task's threads, and 40% of its memory split between them
        def sort_threads = Math.min(4, task.cpus as int)
        def sort_mem = task.memory ? "-m ${((task.memory.toMega() * 0.4) / sort_threads) as int}M" : ''
        """
        minimap2 \\
            -a ${preset} \\
            ${args} \\
            -t ${task.cpus} \\
            --split-prefix ${prefix}.split \\
            ${reference} \\
            ${reads} \\
            | samtools sort \\
                -@ ${sort_threads} \\
                ${sort_mem} \\
                -T ${prefix}.sort_tmp \\
                --write-index \\
                -o ${prefix}.bam##idx##${prefix}.bam.bai \\
                -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        minimap2 \\
            ${preset} \\
            ${args} \\
            -t ${task.cpus} \\
            ${reference} \\
            ${reads} \\
            > ${prefix}.paf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version)
        END_VERSIONS
        """
    }
}
