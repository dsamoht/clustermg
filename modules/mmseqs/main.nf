process MMSEQS_EASYCLUSTER {

    conda "bioconda::mmseqs2=15.6f452"
    if (workflow.containerEngine == 'singularity') {
        container = params.mmseqs_singularity
    } else {
        container = params.mmseqs_docker
    }

    publishDir "${params.outdir}/mmseqs/", mode: 'copy'

    input:
    tuple val(meta), path(aminoacids)
        
    output:
    tuple val(meta), path('cluster_res_cluster.tsv'), emit: clusters_tsv

    script:
    """
    mmseqs easy-cluster ${aminoacids} cluster_res tmp --min-seq-id 0.7 --cov-mode 1 -c 0.9 --threads ${task.cpus}
    """
}
