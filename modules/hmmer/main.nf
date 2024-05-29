process HMMER {

    if (workflow.containerEngine == 'singularity') {
        container = params.hmmer_singularity
    } else {
        container = params.hmmer_docker
    }

    publishDir "${params.outdir}/hmmer", mode: 'copy'

    input:
    path genes
    tuple val(type), path(profile)

    output:
    path 'hmmer_pfam.out', emit: hmmerOutputFilePfam, optional: true
    path 'hmmer_dom-table_pfam.txt', emit: hmmerDomTablePfam, optional: true
    path 'hmmer_kegg.out', emit: hmmerOutputFileKegg, optional: true
    path 'hmmer_dom-table_kegg.txt', emit: hmmerDomTableKegg, optional: true

    script:
    """
    hmmsearch -E 0.001 -o hmmer_${type}.out --domtbl hmmer_dom-table_${type}.txt ${profile} ${genes}
    """
}