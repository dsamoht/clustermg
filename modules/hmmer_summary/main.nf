process HMMER_SUMMARY {

    conda "conda-forge::pandas=2.2.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    publishDir "${params.outdir}/hmmer", mode: 'copy'

    input:
    tuple val(meta), path(hmmerDomTable)
    path koList
    tuple val(meta), path(diamond_result)

    output:
    tuple val(meta), path('genes_annot_summary.tsv'), emit: hmmerSummary

    script:
    hmmerDomTablePfam = hmmerDomTable.grep(~/.*pfam.*/)
    hmmerDomTableKegg = hmmerDomTable.grep(~/.*kegg.*/)
    if(hmmerDomTablePfam.isEmpty()) {
        hmmerDomTablePfam = hmmerDomTablePfam.join('')
    } else {
        hmmerDomTablePfam = '-p ' + hmmerDomTablePfam.join('')
    }
    if(hmmerDomTableKegg.isEmpty()) {
        hmmerDomTableKegg = hmmerDomTableKegg.join('')
    } else {
        hmmerDomTableKegg = '-k ' + hmmerDomTableKegg.join('')
    }
    """
    genes_annot_summary.py ${hmmerDomTablePfam} ${hmmerDomTableKegg} -l ${koList} -d ${diamond_result}
    """
}