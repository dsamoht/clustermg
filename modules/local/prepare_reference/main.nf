process PREPARE_REFERENCE {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(contigs, stageAs: 'contigs/*')
    path(gffs, stageAs: 'gff/*')
    path(contig2mag, stageAs: 'contig2mag/*')
    path(mags_summary)

    output:
    tuple val(meta), path("reference.fa")                         , emit: reference
    tuple val(meta), path("link_query.fa"), path("link_target.fa"), emit: link_sequences
    tuple val(meta), path("contigs.bed")                          , emit: contigs_bed
    tuple val(meta), path("reference_genes.saf")                  , emit: saf
    path "versions.yml"                                           , emit: versions

    script:
    """
    prepare_reference.py \\
        --contigs contigs/* \\
        --gff gff/* \\
        --contig2mag contig2mag/* \\
        --mags-summary ${mags_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
