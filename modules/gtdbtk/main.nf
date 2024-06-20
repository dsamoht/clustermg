process GTDBTK {

    conda "bioconda::gtdbtk=2.3.2"
    if (workflow.containerEngine == 'singularity') {
        container = params.gtdbtk_singularity
    } else {
        container = params.gtdbtk_docker
    }
    
    errorStrategy 'ignore'

    publishDir "${params.outdir}/gtdbtk", mode: 'copy'

    input:
    tuple val(meta), path("bins/*")
    path gtdbtk_db

    output:
    tuple val(meta), path("gtdbtk.*.summary.tsv"), emit: summary

    script:
    """
    export GTDBTK_DATA_PATH=${gtdbtk_db}
    gtdbtk classify_wf --genome_dir bins --out_dir . --skip_ani_screen --extension .fa --cpus ${task.cpus}
    """
}