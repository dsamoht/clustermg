process PREPARE_SAMPLE_DB {

    conda "conda-forge::pandas=2.2.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path('*_gene_info.tsv'), emit: gene_info

    script:
    def name = meta.name
    def fcounts = data[0]
    def faa = data[1]
    """
    build_sample_gene_db.py --name ${name} --faa ${faa} --fcounts ${fcounts}
    """
}