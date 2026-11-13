process MINIMAP2_INDEX {

    label 'process_medium'

    tag "$preset"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.30--h577a1d6_0' :
        'biocontainers/minimap2:2.30--h577a1d6_0' }"

    input:
    tuple val(meta), val(preset), path(reference)

    output:
    tuple val(meta), val(preset), path("*.mmi"), emit: index
    path "versions.yml"                         , emit: versions

    script:
    // An index holds the k-mer settings of its preset, so each preset in use gets its own.
    // A reference larger than the batch size (-I, 8G by default) is indexed in parts, which
    // MINIMAP2_ALIGN merges back with --split-prefix.
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -x ${preset} \\
        ${args} \\
        -t ${task.cpus} \\
        -d reference.${preset}.mmi \\
        ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
