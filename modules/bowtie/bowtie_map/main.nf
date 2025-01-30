process BOWTIE_MAP {

    conda "bioconda::bowtie2=2.4.2"
    if (workflow.containerEngine == 'singularity') {
        container = params.bowtie_singularity
    } else {
        container = params.bowtie_docker
    }

    publishDir "${params.outdir}/bowtie2/", mode: 'copy'

    input:
    tuple val(meta), val(data), path(index)

    output:
    tuple val(meta), val(data), path('*.sam'), emit: sam

    script:
    def name = meta.name
    def input = "-1 \"${data[1]}\" -2 \"${data[2]}\""
    """
    bowtie2 \\
        -p "${task.cpus}" \\
        -x ${name} \\
        $input \\
        > ${name}.sam
    """
}
