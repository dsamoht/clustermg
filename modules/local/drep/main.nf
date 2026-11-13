process DREP {

    label 'process_high'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.7.1--pyhdfd78af_0' :
        'biocontainers/drep:3.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(passing), path(genome_info)
    path(mags, stageAs: 'mags/*')

    output:
    tuple val(meta), path("{Cdb,Ndb,Sdb,Wdb,Widb,genomeInformation}.csv"), emit: tables
    path "versions.yml"                                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Every MAG of the run is staged; dRep is given the ones that passed the filter
    """
    sed 's|^|mags/|' ${passing} > genomes.txt

    dRep dereplicate drep_out \\
        -g genomes.txt \\
        --genomeInfo ${genome_info} \\
        -p ${task.cpus} \\
        ${args}

    # The tables worth keeping, moved out of dRep's work directory so that they are published
    # on their own: clusters, secondary ANI comparisons, scores, winners and genome statistics.
    # Mdb (all-vs-all Mash distances) grows with the square of the number of MAGs, and Bdb only
    # lists the paths dRep read the genomes from.
    mv drep_out/data_tables/{Cdb,Ndb,Sdb,Wdb,Widb,genomeInformation}.csv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep -h 2>&1 | sed -n 's/.*dRep v\\([0-9.]*\\).*/\\1/p' | head -n 1)
        skani: \$(skani --version 2>&1 | sed 's/skani //')
    END_VERSIONS
    """
}
