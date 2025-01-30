#!/usr/bin/env nextflow

include { BOWTIE_BUILD } from '../../modules/bowtie/bowtie_build'
include { BOWTIE_MAP   } from '../../modules/bowtie/bowtie_map'
include { SAMTOOLS     } from '../../modules/bowtie/samtools/'


workflow BOWTIE {

    take:
    sample_input

    main:
    
    sample_input
        .filter { meta, files ->
            files[4] == false
        }
        .set { ch_to_bowtie }

    to_bowtie_map_ch = BOWTIE_BUILD(ch_to_bowtie)
    to_samtools_ch = BOWTIE_MAP(to_bowtie_map_ch)
    ch_aligned_samples = SAMTOOLS(to_samtools_ch)
        .map { meta, files, sorted_bam ->
            return [ meta, [files[0], files[1], files[2], files[3], sorted_bam ] ]
        }

    sample_input
        .filter { meta, files ->
            files[4] != false
        }
        .set { ch_with_bam }

    ch_all_samples = ch_with_bam.mix(ch_aligned_samples)

    emit:
    ch_all_samples
}