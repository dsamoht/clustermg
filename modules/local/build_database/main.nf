process BUILD_DATABASE {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_python-duckdb:d744b89f644c7927"

    input:
    tuple val(meta), path(assemblies)
    path(contig_names, stageAs: 'contig_names/*')
    path(contig2mag, stageAs: 'contig2mag/*')
    path(mags_summary)
    path(gffs, stageAs: 'gff/*')
    path(gene_clusters)
    path(catalog_faa)
    path(catalog_fna)
    path(counts, stageAs: 'counts/*')

    output:
    tuple val(meta), path("clustermg.duckdb")   , emit: database
    tuple val(meta), path("parquet")            , emit: parquet
    tuple val(meta), path("database_checks.tsv"), emit: checks
    path "versions.yml"                         , emit: versions

    script:
    // `counts` is empty without read libraries: the count tables are created empty
    def counts_args = counts ? '--counts counts' : ''
    """
    build_database.py \\
        --assemblies ${assemblies} \\
        --contig-names contig_names/* \\
        --contig2mag contig2mag/* \\
        --mags-summary ${mags_summary} \\
        --gff gff/* \\
        --gene-clusters ${gene_clusters} \\
        --catalog-faa ${catalog_faa} \\
        --catalog-fna ${catalog_fna} \\
        ${counts_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        duckdb: \$(python3 -c 'import duckdb; print(duckdb.__version__)')
    END_VERSIONS
    """
}
