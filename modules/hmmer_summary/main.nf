process HMMER_SUMMARY {

    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    publishDir "${params.outdir}/hmmer", mode: 'copy'

    input:
    path hmmerDomTablePfam
    path hmmerDomTableKegg

    output:
    path 'contig_annotation.tsv', emit: hmmerSummary


    script:
    """
    python hmmer_summary.py -p ${hmmerDomTablePfam} -k ${hmmerDomTableKegg}
    """
}