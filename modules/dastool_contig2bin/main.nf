process DASTOOL_CONTIG2BIN {

    if (workflow.containerEngine == 'singularity') {
        container = params.dastool_singularity
    } else {
        container = params.dastool_docker
    }

    errorStrategy 'ignore'

    publishDir "${params.outdir}/dastool", mode: 'copy'

    input:
    tuple val(meta), path(bins)
    val software

    output:
    tuple val(meta), path("${software}_contigs2bins.tsv"), emit: contigs2bins, optional: true

    script:
    """
    Fasta_to_Contig2Bin.sh -e fa > ${software}_contigs2bins.tsv
    """
}