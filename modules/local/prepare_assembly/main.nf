process PREPARE_ASSEMBLY {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(assembly, stageAs: 'input/*'), path(bins, stageAs: 'input/bins'), path(checkm2, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.contigs.fa")             , emit: contigs
    tuple val(meta), path("mags/*.fa", arity: '0..*') , emit: mags
    tuple val(meta), path("*.contig2mag.tsv")         , emit: contig2mag
    tuple val(meta), path("*.mag_info.tsv")           , emit: mag_info
    tuple val(meta), path("*.contig_names.tsv")       , emit: contig_names
    path "versions.yml"                               , emit: versions

    script:
    // The assembly name is not a file prefix: it becomes part of every contig, gene and MAG
    // id, so it always is `meta.id`.
    def bins_args = bins ? "--bins ${bins} --checkm2 ${checkm2}" : ''
    """
    prepare_assembly.py \\
        --name ${meta.id} \\
        --assembly ${assembly} \\
        ${bins_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
