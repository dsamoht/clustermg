process FEATURECOUNTS {

    conda "bioconda::subread=2.0.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.subread_singularity
    } else {
        container = params.subread_docker
    }

    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    tuple val(meta), val(data)
    val read_type

    output:
    tuple val(meta), path("*featureCounts.txt"), emit: counts

    script:
    if ("${read_type}" == "long") {
        options = "-L -t CDS,gene -g ID -s 0"
    } else {
        options = "-p -t CDS,gene -g ID -s 0"
    }
    def genes_gff = data[0]
    def sorted_bam = data[1]
    def name = meta.name
    """
    featureCounts ${options} -a ${genes_gff} -o ${name}_featureCounts.txt ${sorted_bam}
    """
}
