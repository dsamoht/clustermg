process PREPARE_DB_DB {

    conda "conda-forge::pandas=2.2.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path('*_gene_info.tsv'), emit: gene_info

    script:
    def name = meta.name
    def faa = data[0]
    """
    build_db_gene_db.py --name ${name} --faa ${faa}
    """
}