process PARTITION_GENES {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(faa)
    path(final_contig2mag)

    output:
    tuple val(meta), path("*.mag_genes.faa")  , emit: mag_genes
    tuple val(meta), path("*.other_genes.faa"), emit: other_genes
    path "versions.yml"                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    partition_genes.py \\
        --faa ${faa} \\
        --final-contig2mag ${final_contig2mag} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
