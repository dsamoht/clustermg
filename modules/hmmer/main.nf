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
    if(params.profilePfam != '') {
    path 'hmmer_pfam.out', emit: hmmerOutputFilePfam
    path 'hmmer_dom-table_pfam.txt', emit: hmmerDomTablePfam
    } else {
        val '', emit: hmmerDomTablePfam
    }
    if(params.profileHmmKegg != '') {
    path 'hmmer_kegg.out', emit: hmmerOutputFileKegg
    path 'hmmer_dom-table_kegg.txt', emit: hmmerDomTableKegg
    } else {
        val '', emit: hmmerDomTableKegg
    }

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