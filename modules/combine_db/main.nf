process COMBINE_DB {

    conda "conda-forge::pandas=2.2.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    publishDir "${params.outdir}/gene_db", mode: 'copy'

    input:
    path db_files

    output:
    path('gene_db.tsv'), emit: combined_sample_db

    script:
    """
    combine_dbs.py --dbs ${db_files}
    """
}
