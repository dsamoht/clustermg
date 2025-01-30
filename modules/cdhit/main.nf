process CDHIT {

    conda "bioconda::cd-hit=4.8.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.cdhit_singularity
    } else {
        container = params.cdhit_docker
    }

    publishDir "${params.outdir}/cdhit", mode: 'copy'

    input:
    tuple val(meta), path(aminoacids)
        
    output:
    tuple val(meta), path('cdhit_c70.clstr'), emit: clusters_tsv

    script:
    """
    cd-hit -i ${aminoacids} -o cdhit_c70 -c 0.70 -aS 0.90 -g 1 -G 0 -d 0 -T ${task.cpus} -M 0
    """
}
