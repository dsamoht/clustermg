process SAMTOOLS_BEDCOV {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(regions)

    output:
    tuple val(meta), path("*.bedcov.tsv"), emit: bedcov
    path "versions.yml"                  , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // DNA is counted on every contig and on the stretches linked to a MAG, in one pass: aligned
    // M/=/X bases (-j), bases at depth 1 or more (-d 1) and reads (-c). bedcov runs without
    // --max-depth, which leaves the depth uncapped (0 would cap it). Its default filter leaves
    // unmapped, secondary, QC-failed and duplicate alignments out; -G 0x800 adds supplementary
    // ones, so that a chimeric long read is counted once, on its primary alignment, as RNA is.
    // No mapping quality threshold applies: a read placed at random between identical copies is
    // still counted once.
    """
    samtools bedcov -c -j -d 1 -G 0x800 ${regions} ${bam} > ${prefix}.bedcov.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
