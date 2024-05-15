process HMMER {

    if (workflow.containerEngine == 'singularity') {
        container = params.hmmer_singularity
    } else {
        container = params.hmmer_docker
    }

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path genes
    path profileHmm

    output:
    path 'hmmer.out', emit: hmmerOutputFile
    path 'hmmer_dom-table.txt', emit: hmmerDomainTable

    script:
    """
    hmmsearch -E 0.001 -o hmmer.out --domtbl hmmer_dom-table.txt ${profileHmm} ${genes}
    """
}