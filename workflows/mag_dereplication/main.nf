#!/usr/bin/env nextflow

include { DREP           } from '../../modules/local/drep'
include { FILTER_MAGS    } from '../../modules/local/filter_mags'
include { SUMMARIZE_MAGS } from '../../modules/local/summarize_mags'

workflow MAG_DEREPLICATION {

    take:
    ch_mags       // channel: [ meta, [ MAG fasta ] ], one entry per assembly, an empty list without bins
    ch_mag_info   // channel: [ meta, mag_info.tsv ], one entry per assembly
    ch_contig2mag // channel: [ meta, contig2mag.tsv ], one entry per assembly

    main:
    ch_versions = channel.empty()

    // The fan-ins below are sorted: the order the assemblies finish in changes from one
    // run to the next, and it would change the task hashes with it, so `-resume` would
    // redo every run-wide step. toSortedList() also emits an empty list when there is
    // nothing to collect, so that a run without a single MAG still goes through.
    ch_all_mags = ch_mags
        .flatMap { _meta, mags -> mags }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_all_mag_info = ch_mag_info
        .map { _meta, mag_info -> mag_info }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_all_contig2mag = ch_contig2mag
        .map { _meta, contig2mag -> contig2mag }
        .toSortedList { a, b -> a.name <=> b.name }

    //
    // MODULE: keep the MAGs good enough to be dereplicated
    //
    FILTER_MAGS(ch_all_mag_info.map { mag_info -> [ [ id: 'mags' ], mag_info ] })
    ch_versions = ch_versions.mix(FILTER_MAGS.out.versions)

    //
    // MODULE: dereplicate them
    //
    // dRep compares genomes: with fewer than two MAGs past the filter there is nothing to
    // compare, the one that passed (if any) is final on its own, and SUMMARIZE_MAGS gets
    // no dRep table.
    ch_to_dereplicate = FILTER_MAGS.out.passing
        .join(FILTER_MAGS.out.genome_info)
        .filter { _meta, passing, _genome_info -> passing.readLines().size() > 1 }

    DREP(ch_to_dereplicate, ch_all_mags)
    ch_versions = ch_versions.mix(DREP.out.versions)

    //
    // MODULE: give every MAG its status and collect the final ones
    //
    SUMMARIZE_MAGS(
        FILTER_MAGS.out.quality,
        DREP.out.tables.map { _meta, tables -> tables }.ifEmpty([]),
        ch_all_contig2mag,
        ch_all_mags
    )
    ch_versions = ch_versions.mix(SUMMARIZE_MAGS.out.versions)

    emit:
    summary          = SUMMARIZE_MAGS.out.summary          // channel: [ meta, mags_summary.tsv ]
    final_contig2mag = SUMMARIZE_MAGS.out.final_contig2mag // channel: [ meta, final_contig2mag.tsv ]
    final_mags       = SUMMARIZE_MAGS.out.final_mags       // channel: [ meta, [ final MAG fasta ] ]
    drep_tables      = DREP.out.tables                     // channel: [ meta, [ dRep data tables ] ]
    versions         = ch_versions                         // channel: [ versions.yml ]
}
