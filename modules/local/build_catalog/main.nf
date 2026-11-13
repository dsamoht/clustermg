process BUILD_CATALOG {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(gffs, stageAs: 'predicted/*'), path(fnas, stageAs: 'predicted/*'), path(faas, stageAs: 'predicted/*')
    path(mag_clusters, stageAs: 'clusters/mag/*')
    path(assigned, stageAs: 'clusters/assigned/*')
    path(other_clusters, stageAs: 'clusters/other/*')
    path(contig2mag, stageAs: 'contig2mag/*')
    path(mags_summary)

    output:
    tuple val(meta), path("gene_catalog.faa")      , emit: faa
    tuple val(meta), path("gene_catalog.fna")      , emit: fna
    tuple val(meta), path("gene_catalog.gff")      , emit: gff
    tuple val(meta), path("gene_clusters.tsv")     , emit: clusters
    tuple val(meta), path("gene_catalog_stats.tsv"), emit: stats
    path "versions.yml"                            , emit: versions

    script:
    // The two clustering steps are skipped when they have no sequence to cluster: without
    // a final MAG, and when every other gene joined a final-MAG representative
    def mag_args   = mag_clusters ? "--mag-clusters ${mag_clusters}" : ''
    def other_args = other_clusters ? "--other-clusters ${other_clusters}" : ''
    """
    build_catalog.py \\
        ${mag_args} \\
        --assigned ${assigned} \\
        ${other_args} \\
        --faa predicted/*.faa \\
        --fna predicted/*.fna \\
        --gff predicted/*.gff \\
        --contig2mag contig2mag/* \\
        --mags-summary ${mags_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
