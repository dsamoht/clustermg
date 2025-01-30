process MAKE_MATRIX {

    conda "conda-forge::pandas=2.2.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.pandas_singularity
    } else {
        container = params.pandas_docker
    }

    publishDir "${params.outdir}/gene_db", mode: 'copy'

    input:
    path clstr
    path gene_db

    output:
    path('mmseqs_cluster_rpkm.tsv'), emit: mmseqs_cluster_rpkm

    script:
    """
    build_matrix_from_clustering.py --clstr ${clstr} --genedb ${gene_db}
    """
}
