process SORTMERNA {

    label 'process_medium'

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.7--hdbdd923_1' :
        'biocontainers/sortmerna:4.3.7--hdbdd923_1' }"

    input:
    tuple val(meta), path(reads)
    path(fastas, stageAs: 'rrna/*')
    path(index)

    output:
    tuple val(meta), path("*.non_rRNA.fastq.gz"), emit: reads
    tuple val(meta), path("*.sortmerna.log")    , emit: log
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // a single database comes as one path rather than a list
    def refs = (fastas instanceof List ? fastas : [ fastas ]).collect { fasta -> "--ref ${fasta}" }.join(' ')
    // Both mates are removed when either is rRNA (--paired_in), and the mates kept go to two
    // files (--out2). The rRNA reads themselves are not kept.
    def mates = meta.single_end ? '' : "--reads ${reads[1]} --paired_in --out2"
    def rename = meta.single_end ?
        "mv non_rRNA_reads.f*q.gz ${prefix}.non_rRNA.fastq.gz" :
        "mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz && mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz"
    """
    sortmerna \\
        ${refs} \\
        --reads ${reads[0]} \\
        ${mates} \\
        --threads ${task.cpus} \\
        --workdir . \\
        --idx-dir ${index} \\
        --index 0 \\
        --aligned rRNA_reads \\
        --fastx \\
        --other non_rRNA_reads \\
        ${args}

    ${rename}
    mv rRNA_reads.log ${prefix}.sortmerna.log
    rm -f rRNA_reads*.f*q.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS
    """
}
