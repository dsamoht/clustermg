process SUMMARIZE_MAGS {

    label 'process_single'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "docker.io/dsamoht/bio-utils@sha256:f0cad0d32d8d8fac7bb971736f158200cf19b4817dd796bd9d76240a054bacf2"

    input:
    tuple val(meta), path(quality)
    path(drep_tables, stageAs: 'drep/*')
    path(contig2mag, stageAs: 'contig2mag/*')
    path(mags, stageAs: 'mags/*')

    output:
    tuple val(meta), path("mags_summary.tsv")                , emit: summary
    tuple val(meta), path("final_contig2mag.tsv")            , emit: final_contig2mag
    tuple val(meta), path("*.fa", arity: '0..*')             , emit: final_mags
    path "versions.yml"                                      , emit: versions

    script:
    // `drep_tables` is empty when dRep did not run (fewer than two MAGs passed the filter),
    // and so is `mags` when no assembly came with bins
    def drep_args = drep_tables ? "--cdb drep/Cdb.csv --sdb drep/Sdb.csv --wdb drep/Wdb.csv" : ''
    def mags_args = mags ? "--mags mags/*" : ''
    """
    summarize_mags.py \\
        --quality ${quality} \\
        --contig2mag contig2mag/* \\
        ${drep_args} \\
        ${mags_args} \\
        --final-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
