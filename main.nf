#!/usr/bin/env nextflow

include { BOWTIE               } from './workflows/bowtie'
include { CHECK_FORMAT         } from './workflows/check_format'
include { COMBINE_DB           } from './modules/combine_db'
include { DISPATCH             } from './workflows/dispatch'
include { FEATURECOUNTS        } from './modules/featurecounts'
include { MAKE_MATRIX          } from './modules/make_matrix'
include { MMSEQS_EASYCLUSTER   } from './modules/mmseqs'
include { PREPARE_DB_DB        } from './modules/prepare_db/prepare_db_db'
include { PREPARE_SAMPLE_DB    } from './modules/prepare_db/prepare_sample_db'
include { PRODIGAL             } from './modules/prodigal'
include { SPLIT_TO_FASTA       } from './modules/split_to_fasta'


info = """
         __              __                               
  _____ / /__  __ _____ / /_ ___   _____ ____ ___   ____ _
 / ___// // / / // ___// __// _ \\ / ___// __ `__ \\ / __ `/
/ /__ / // /_/ /(__  )/ /_ /  __// /   / / / / / // /_/ / 
\\___//_/ \\__,_//____/ \\__/ \\___//_/   /_/ /_/ /_/ \\__, /  
                                                 /____/
     
     Github: https://github.com/dsamoht/clustermg
     Version: v1.0.0

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

if( params.help ) {

log.info info
    exit 0
}

log.info info

workflow {

    (ch_sample, ch_db) = DISPATCH() | CHECK_FORMAT

    ch_prodigal_input = ch_sample
        .map{ meta, files -> return [ meta, files[0] ] }

    ch_prodigal_output = PRODIGAL(ch_prodigal_input)
    
    ch_bam_for_featurecounts = BOWTIE(ch_sample)
        .map { meta, files -> 
          return [ meta, files[3] ] }
   
    ch_featurecounts_in = ch_prodigal_output.genes_gff
        .concat(ch_bam_for_featurecounts)
        .groupTuple()
       
    ch_featurecounts_out = FEATURECOUNTS(ch_featurecounts_in, "short_reads")

    ch_to_construct_gene_db = ch_featurecounts_out.counts
        .concat(ch_prodigal_output.genes_faa)
        .groupTuple()

    ch_formatted_sample_db = PREPARE_SAMPLE_DB(ch_to_construct_gene_db)
        .map{ meta, files -> return files }
        .collect()

    ch_formatted_db_db = PREPARE_DB_DB(ch_db)
        .map{ meta, files -> return files }
        .collect()

    ch_formatted_dbs = ch_formatted_sample_db
        .concat(ch_formatted_db_db)
        .collect()

    ch_combined_database = COMBINE_DB(ch_formatted_dbs)

    (ch_sample_db, ch_db_db) = SPLIT_TO_FASTA(ch_combined_database)

    ch_sample_db = ch_sample_db
        .map { db -> 
                def meta = [:]
                meta.name = 'sample'
                return [meta, db]
             }

    ch_sample_mmseqs = MMSEQS_EASYCLUSTER(ch_sample_db)
        .map { meta, file -> return file }

    MAKE_MATRIX(ch_sample_mmseqs, ch_combined_database)

}
