#!/usr/bin/env nextflow

include { SEQKIT                   } from '../../modules/seqkit/'


workflow CHECK_FORMAT {

    take:
    sample_input
    db_input

    main:
    ch_input = sample_input.mix(db_input)
        .map{meta, files -> return [meta, files[0]]}
    
    SEQKIT(ch_input).seq_type
                .map { meta, seq_type ->
                if (seq_type.text.trim() == "Protein" && meta.type == "sample") {
                    exit 1, "Error: '$meta.name' fasta file contains protein sequences"
                    }
                if (seq_type.text.trim() == "DNA" && meta.type == "db") {
                    exit 1, "Error: '$meta.name' fasta file contains DNA sequences"
                    }
                }

    emit:
    sample_input
    db_input

}
