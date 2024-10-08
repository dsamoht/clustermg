process FEATURECOUNTS {

    conda "bioconda::subread=2.0.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.subread_singularity
    } else {
        container = params.subread_docker
    }

    errorStrategy 'ignore'
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    tuple val(meta), path(genes_gff)
    tuple val(meta), path(sorted_bam)
    val read_type

    output:
    tuple val(meta), path("*featureCounts.txt"), emit: counts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary

    script:
    if ("${read_type}" == "long") {
        options = "-L -t CDS,gene -g ID -s 0"
    } else {
        options = "-p -t CDS,gene -g ID -s 0"
    }
    def genes = genes_gff.join(' ')
    """
    cat ${genes} > global.gff 
    featureCounts ${options} -a global.gff -o featureCounts.txt ${sorted_bam}
    """
}
