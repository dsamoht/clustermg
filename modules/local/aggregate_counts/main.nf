process AGGREGATE_COUNTS {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(library_sheet)
    path(flagstats, stageAs: 'flagstat/*')
    path(bedcovs, stageAs: 'bedcov/*')
    path(featurecounts, stageAs: 'featurecounts/*')
    path(sortmerna_logs, stageAs: 'sortmerna/*')
    path(coverm, stageAs: 'coverm/*')
    path(contig_links)
    path(contig2mag, stageAs: 'contig2mag/*')
    path(mags_summary)
    path(gene_clusters)
    path(saf)

    output:
    tuple val(meta), path("libraries.tsv")        , emit: libraries
    tuple val(meta), path("contig_dna_counts.tsv"), emit: contig_dna_counts
    tuple val(meta), path("link_dna_counts.tsv")  , emit: link_dna_counts
    tuple val(meta), path("mag_dna_profile.tsv")  , emit: mag_dna_profile
    tuple val(meta), path("final_mag_dna.tsv")    , emit: final_mag_dna
    tuple val(meta), path("fine_group_dna.tsv")   , emit: fine_group_dna
    tuple val(meta), path("gene_rna_counts.tsv")  , emit: gene_rna_counts
    tuple val(meta), path("gene_expression.tsv")  , emit: gene_expression
    path "versions.yml"                           , emit: versions

    script:
    // DNA libraries come with bedcov tables and RNA libraries with featureCounts tables, so
    // either list is empty when the run has no library of that molecule; so is `coverm`
    // without a MAG
    def bedcov_args = bedcovs ? '--bedcov bedcov/*' : ''
    def featurecounts_args = featurecounts ? '--featurecounts featurecounts/*.featureCounts.tsv' : ''
    def coverm_args = coverm ? "--coverm ${coverm}" : ''
    // and without rRNA removal, there is no SortMeRNA log
    def sortmerna_args = sortmerna_logs ? '--sortmerna-log sortmerna/*' : ''
    """
    aggregate_counts.py \\
        --libraries ${library_sheet} \\
        --flagstat flagstat/* \\
        ${bedcov_args} \\
        ${featurecounts_args} \\
        ${coverm_args} \\
        ${sortmerna_args} \\
        --contig-links ${contig_links} \\
        --contig2mag contig2mag/* \\
        --mags-summary ${mags_summary} \\
        --gene-clusters ${gene_clusters} \\
        --saf ${saf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
