#!/usr/bin/env nextflow

workflow DISPATCH {

    ch_input_rows = Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() == 7) {
                    def name = row.name
                    def type = row.type
                    def fasta = row.fasta ? file(row.fasta, checkIfExists: true) : false
                    def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : false
                    def sr2 = row.short_reads_2 ? file(row.short_reads_2, checkIfExists: true) : false
                    def lr = row.long_reads ? file(row.long_reads, checkIfExists: true) : false
                    def bam = row.bam ? file(row.bam, checkIfExists: true) : false
                    return [ name, type, fasta, sr1, sr2, lr, bam ]
     
                } else {
                    exit 1, "Error in ${params.input}. Each row must contain 7 columns."
                }
        }
    
    ch_sample_input = ch_input_rows
        .filter { it[1] == 'sample' }
        .map { name, type, fasta, sr1, sr2, lr, bam -> 
                def meta = [:]
                meta.name = name
                meta.type = type
                return [ meta, [ fasta, sr1, sr2, lr, bam ] ]
              }
    
    ch_db_input = ch_input_rows
        .filter { it[1] == 'db' }
        .map { name, type, fasta, sr1, sr2, lr, bam -> 
                def meta = [:]
                meta.name = name
                meta.type = type
                return  [ meta, [ fasta ] ]
              }

    emit:
    sample_input = ch_sample_input
    db_input = ch_db_input

}
