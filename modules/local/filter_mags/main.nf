process FILTER_MAGS {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(mag_info, stageAs: 'mag_info/*')

    output:
    tuple val(meta), path("mag_quality.tsv") , emit: quality
    tuple val(meta), path("genome_info.csv") , emit: genome_info
    tuple val(meta), path("passing_mags.txt"), emit: passing
    path "versions.yml"                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    filter_mags.py ${args} mag_info/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
