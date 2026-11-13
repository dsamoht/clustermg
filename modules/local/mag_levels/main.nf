process MAG_LEVELS {

    label 'process_single'

    tag "$meta.id"

    // dRep's own container: the tree is rebuilt with the scipy dRep clusters with
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.7.1--pyhdfd78af_0' :
        'biocontainers/drep:3.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(drep_tables, stageAs: 'drep/*')

    output:
    tuple val(meta), path("mag_levels.tsv"), emit: levels
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mag_levels.py \\
        --ndb drep/Ndb.csv \\
        --cdb drep/Cdb.csv \\
        --sdb drep/Sdb.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
    END_VERSIONS
    """
}
