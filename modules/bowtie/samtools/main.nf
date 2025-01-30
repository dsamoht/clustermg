process SAMTOOLS {

    conda "bioconda::samtools=1.18"
    if (workflow.containerEngine == 'singularity') {
        container = params.samtools_singularity
    } else {
        container = params.samtools_docker
    }

    publishDir "${params.outdir}/samtools", mode: 'copy'

    input:
    tuple val(meta), val(data), path(sam)

    output:
    tuple val(meta), val(data), path('*.sorted.bam'), emit: sorted_bam

    script:
    def name = meta.name
    """
    samtools view -bS ${sam} | samtools sort -o ${name}.sorted.bam -
    """
}
