process ASSIGN_GENES {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(proteins, stageAs: 'proteins/*')
    path(hits)

    output:
    tuple val(meta), path("assigned.tsv")  , emit: assigned
    tuple val(meta), path("unassigned.faa"), emit: unassigned
    path "versions.yml"                    , emit: versions

    script:
    // `hits` is empty when there was no final-MAG representative to search against, and
    // `proteins` when every gene of the run sits on a final-MAG contig
    def hits_args = hits ? "--hits ${hits}" : ''
    def proteins_args = proteins ? "proteins/*" : ''
    """
    assign_genes.py ${hits_args} ${proteins_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
