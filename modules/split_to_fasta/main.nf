process SPLIT_TO_FASTA {

    conda "conda-forge::pandas=2.2.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    publishDir "${params.outdir}/gene_db", mode: 'copy'

    input:
    path(combined_db)

    output:
    path('sample.faa'), emit: sample_fasta
    path('db.faa'), emit: db_fasta

    script:
    """
    gene_db_to_fasta.py --db ${combined_db}
    """
}