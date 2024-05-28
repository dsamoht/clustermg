process HMMER {

    if (workflow.containerEngine == 'singularity') {
        container = params.hmmer_singularity
    } else {
        container = params.hmmer_docker
    }

    publishDir "${params.outdir}/hmmer", mode: 'copy'

    input:
    path genes
    path profileHmmPfam
    path profileHmmKegg

    output:
    path 'hmmer_*.out', emit: hmmerOutputFile
    path 'hmmer_dom-table_*.txt', emit: hmmerDomTable

    script:
    """
    if [[ !  -z ${profileHmmPfam} ]]
    then
    hmmsearch -E 0.001 -o hmmer_pfam.out --domtbl hmmer_dom-table_pfam.txt ${profileHmmPfam} ${genes}
    fi
    if [[ !  -z ${profileHmmKegg} ]]
    then
    hmmsearch -E 0.001 -o hmmer_kegg.out --domtbl hmmer_dom-table_kegg.txt ${profileHmmKegg} ${genes}
    fi
    """
}