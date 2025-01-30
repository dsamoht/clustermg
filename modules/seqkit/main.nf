process SEQKIT {

    conda "bioconda::seqkit=2.5.1"
    if (workflow.containerEngine == 'singularity') {
        container = params.seqkit_singularity
    } else {
        container = params.seqkit_docker
    }

    publishDir "${params.outdir}/seqkit/", mode: 'copy', pattern: '*_fasta_stats.tsv'

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path('*_fasta_stats.tsv')
    tuple val(meta), path('*_seq_type.txt'), emit: seq_type

    script:
    def fasta = data[0]
    def name = meta.name
    """
    seqkit stats -aTb ${fasta} > ${name}_fasta_stats.tsv
    seq_type=\$(tail -1 ${name}_fasta_stats.tsv | awk '{print \$3}')
    echo \$seq_type > ${name}_seq_type.txt
    """
}
