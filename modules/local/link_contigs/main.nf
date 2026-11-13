process LINK_CONTIGS {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(contigs_bed), path(paf)
    path(contig2mag, stageAs: 'contig2mag/*')

    output:
    tuple val(meta), path("contig_links.tsv"), emit: links
    tuple val(meta), path("regions.bed")     , emit: regions
    path "versions.yml"                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    // `paf` is empty without a passing MAG: nothing to link to, the regions are the contigs
    def paf_args = paf ? "--paf ${paf}" : ''
    """
    link_contigs.py \\
        ${paf_args} \\
        --contigs-bed ${contigs_bed} \\
        --contig2mag contig2mag/* \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
