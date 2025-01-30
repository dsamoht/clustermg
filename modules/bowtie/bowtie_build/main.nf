process BOWTIE_BUILD {

    conda "bioconda::bowtie2=2.4.2"
    if (workflow.containerEngine == 'singularity') {
        container = params.bowtie_singularity
    } else {
        container = params.bowtie_docker
    }

    publishDir "${params.outdir}/bowtie2/", mode: 'copy'

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path(data), path('*.bt2'), emit: fasta_index

    script:
    def fasta = data[0]
    def name = meta.name
    """
    bowtie2-build --threads ${task.cpus} ${fasta} ${name}
    """
}
