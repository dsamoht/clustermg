process PRODIGAL {

    conda "bioconda::prodigal=2.6.3"
    if (workflow.containerEngine == 'singularity') {
        container = params.prodigal_singularity
    } else {
        container = params.prodigal_docker
    }

    publishDir "${params.outdir}/prodigal/", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*genes_cds.gff'), emit: genes_gff
    tuple val(meta), path('*genes.faa'), emit: genes_faa
    tuple val(meta), path('*genes.fna'), emit: genes_fna

    script:
    def name = meta.name
    """
    prodigal -i ${fasta} -f gff -o ${name}.genes.gff -a ${name}.genes.faa -d ${name}.genes.fna -p meta
    grep "CDS" ${name}.genes.gff > ${name}.genes_cds.gff
    """
}