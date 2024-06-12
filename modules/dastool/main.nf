process DASTOOL {

    if (workflow.containerEngine == 'singularity') {
        container = params.dastool_singularity
    } else {
        container = params.dastool_docker
    }

    errorStrategy 'ignore'

    publishDir "${params.outdir}/dastool", mode: 'copy', saveAs: { filename -> new File(filename).getName() }

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(contig2bin)

    output:
    tuple val(meta), path("das-bin*/*bin*.fa"), emit: dasBins, optional: true

    script:
    def contig2binList = contig2bin.join(",")
    def label = contig2bin instanceof List ? "metabat,maxbin" : contig2bin.toString() - "_contigs2bins.tsv"
    """
    DAS_Tool -i ${contig2binList} -l ${label} -c ${assembly} -o das-bin --write_bins --score_threshold=-9999
    """
}