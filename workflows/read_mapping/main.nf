#!/usr/bin/env nextflow

include { AGGREGATE_COUNTS                         } from '../../modules/local/aggregate_counts'
include { BWAMEM2_INDEX                            } from '../../modules/local/bwamem2/index'
include { BWAMEM2_MEM                              } from '../../modules/local/bwamem2/mem'
include { COVERM_GENOME                            } from '../../modules/local/coverm/genome'
include { LINK_CONTIGS                             } from '../../modules/local/link_contigs'
include { MINIMAP2_ALIGN                           } from '../../modules/local/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_CONTIGS } from '../../modules/local/minimap2/align'
include { MINIMAP2_INDEX                           } from '../../modules/local/minimap2/index'
include { PREPARE_REFERENCE                        } from '../../modules/local/prepare_reference'
include { SAMTOOLS_BEDCOV                          } from '../../modules/local/samtools/bedcov'
include { SAMTOOLS_FLAGSTAT                        } from '../../modules/local/samtools/flagstat'
include { SORTMERNA                                } from '../../modules/local/sortmerna'
include { SORTMERNA_INDEX                          } from '../../modules/local/sortmerna/index'
include { SUBREAD_FEATURECOUNTS                    } from '../../modules/local/subread/featurecounts'

//
// Every read library is mapped once, competitively, onto the contigs of all the assemblies.
// A read lands on the copy it matches best, the strain's own where several were assembled, and
// on one of them at random where the copies are identical. Each read sits on one copy only, so
// the copies of a genome or of a gene add up without counting it twice:
//
//   - DNA is counted on every contig, and on the stretches of unbinned or low-quality contigs
//     that align to a passing MAG, which count as part of it;
//   - RNA, its rRNA reads removed first, is counted on every gene, and each gene's count goes to
//     its catalog representative.
//
// Illumina reads are mapped with bwa-mem2, long reads with minimap2.
//
workflow READ_MAPPING {

    take:
    ch_reads          // channel: [ meta, [ fastq ] ], one entry per read library
    ch_contigs        // channel: [ meta, contigs.fa ], one entry per assembly
    ch_gff            // channel: [ meta, gff ], one entry per assembly
    ch_contig2mag     // channel: [ meta, contig2mag.tsv ], one entry per assembly
    ch_mags           // channel: [ meta, [ MAG fasta ] ], one entry per assembly
    ch_mags_summary   // channel: [ meta, mags_summary.tsv ]
    ch_gene_clusters  // channel: [ meta, gene_clusters.tsv ]
    ch_rrna_databases // channel: [ rRNA FASTA ], value, empty when the rRNA is not removed

    main:
    ch_versions = channel.empty()

    // Sorted fan-ins, so that the run-wide tasks keep their hash, and -resume its cache,
    // whatever order the assemblies finish in. Nothing here runs without a read library.
    ch_has_reads = ch_reads.count().filter { n -> n > 0 }

    ch_all_contigs = ch_contigs
        .map { _meta, contigs -> contigs }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_all_gff = ch_gff
        .map { _meta, gff -> gff }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_all_contig2mag = ch_contig2mag
        .map { _meta, contig2mag -> contig2mag }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_all_mags = ch_mags
        .flatMap { _meta, mags -> mags }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_summary = ch_mags_summary.map { _meta, summary -> summary }.first()

    //
    // MODULE: the reference is every contig of every assembly
    //
    PREPARE_REFERENCE(
        ch_has_reads
            .combine(ch_all_contigs.map { contigs -> [ contigs ] })
            .map { _n, contigs -> [ [ id: 'reference' ], contigs ] },
        ch_all_gff,
        ch_all_contig2mag,
        ch_summary
    )
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // MODULES: tie the stretches of unbinned and low-quality contigs to the MAG they copy.
    // Without a MAG that passed the filter there is nothing to align to.
    //
    MINIMAP2_ALIGN_CONTIGS(
        PREPARE_REFERENCE.out.link_sequences
            .filter { _meta, query, target -> query.size() > 0 && target.size() > 0 },
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_CONTIGS.out.versions)

    LINK_CONTIGS(
        PREPARE_REFERENCE.out.contigs_bed
            .join(MINIMAP2_ALIGN_CONTIGS.out.paf, remainder: true)
            .map { meta, bed, paf -> [ meta, bed, paf ?: [] ] },
        ch_all_contig2mag
    )
    ch_versions = ch_versions.mix(LINK_CONTIGS.out.versions)

    //
    // MODULES: remove the rRNA reads of the RNA libraries, the index built once for all of them
    //
    ch_dna = ch_reads.filter { meta, _reads -> meta.molecule == 'dna' }
    ch_rna = ch_reads.filter { meta, _reads -> meta.molecule == 'rna' }

    if (params.skip_rrna_removal) {
        ch_rna_to_map = ch_rna
        ch_rrna_logs = channel.empty()
    }
    else {
        SORTMERNA_INDEX(
            ch_rna.count()
                .filter { n -> n > 0 }
                .combine(ch_rrna_databases.map { databases -> [ databases ] })
                .map { _n, databases -> databases }
        )
        ch_versions = ch_versions.mix(SORTMERNA_INDEX.out.versions)

        SORTMERNA(ch_rna, ch_rrna_databases, SORTMERNA_INDEX.out.index.first())
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())

        ch_rna_to_map = SORTMERNA.out.reads
        ch_rrna_logs = SORTMERNA.out.log
    }

    ch_to_map = ch_dna.mix(ch_rna_to_map)

    //
    // MODULES: long reads, onto a minimap2 index per preset in use. PacBio reads are taken as HiFi.
    //
    def presets = [ ont: 'map-ont', pacbio: 'map-hifi' ]
    ch_long = ch_to_map
        .filter { meta, _reads -> meta.platform != 'illumina' }
        .map { meta, reads -> [ meta + [ preset: presets[meta.platform] ], reads ] }

    MINIMAP2_INDEX(
        ch_long
            .map { meta, _reads -> meta.preset }
            .unique()
            .combine(PREPARE_REFERENCE.out.reference)
            .map { preset, meta, reference -> [ meta, preset, reference ] }
    )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())

    MINIMAP2_ALIGN(
        ch_long
            .map { meta, reads -> [ meta.preset, meta, reads ] }
            .combine(MINIMAP2_INDEX.out.index.map { _meta, preset, index -> [ preset, index ] }, by: 0)
            .map { _preset, meta, reads, index -> [ meta, reads, index ] },
        true
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //
    // MODULES: Illumina reads, DNA and RNA alike, onto one bwa-mem2 index
    //
    ch_short = ch_to_map.filter { meta, _reads -> meta.platform == 'illumina' }

    BWAMEM2_INDEX(
        ch_short.count()
            .filter { n -> n > 0 }
            .combine(PREPARE_REFERENCE.out.reference)
            .map { _n, meta, reference -> [ meta, reference ] }
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    BWAMEM2_MEM(ch_short.combine(BWAMEM2_INDEX.out.index.map { _meta, index -> index }))
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    ch_bams = MINIMAP2_ALIGN.out.bam.mix(BWAMEM2_MEM.out.bam)

    //
    // MODULES: the reads of every library, and the DNA on every contig and linked stretch
    //
    SAMTOOLS_FLAGSTAT(ch_bams)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_BEDCOV(
        ch_bams
            .filter { meta, _bam, _bai -> meta.molecule == 'dna' }
            .combine(LINK_CONTIGS.out.regions.map { _meta, regions -> regions })
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BEDCOV.out.versions.first())

    //
    // MODULE: count the RNA on every gene, strand-aware
    //
    SUBREAD_FEATURECOUNTS(
        ch_bams.filter { meta, _bam, _bai -> meta.molecule == 'rna' },
        PREPARE_REFERENCE.out.saf.first()
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // MODULE: profile every MAG, on its own contigs, in every DNA library
    //
    ch_dna_bams = ch_bams
        .filter { meta, _bam, _bai -> meta.molecule == 'dna' }
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .filter { libraries -> libraries }
        .map { libraries -> [ [ id: 'mags' ], libraries.collect { lib -> lib[1] }, libraries.collect { lib -> lib[2] } ] }

    COVERM_GENOME(ch_dna_bams, ch_all_mags.filter { mags -> mags })
    ch_versions = ch_versions.mix(COVERM_GENOME.out.versions)

    //
    // MODULE: gather every count, and add them up into final MAGs, finer groups and catalog genes
    //
    ch_library_sheet = ch_reads
        .map { meta, _reads ->
            [ meta.id, meta.sample, meta.molecule, meta.platform, meta.single_end ? 'single' : 'paired', meta.strandedness ?: '' ].join('\t')
        }
        .collectFile(name: 'library_sheet.tsv', seed: 'library\tsample\tmolecule\tplatform\tlayout\tstrandedness', newLine: true, sort: true)
        .map { sheet -> [ [ id: 'counts' ], sheet ] }

    AGGREGATE_COUNTS(
        ch_library_sheet,
        SAMTOOLS_FLAGSTAT.out.flagstat.map { _meta, flagstat -> flagstat }.toSortedList { a, b -> a.name <=> b.name },
        SAMTOOLS_BEDCOV.out.bedcov.map { _meta, bedcov -> bedcov }.toSortedList { a, b -> a.name <=> b.name },
        SUBREAD_FEATURECOUNTS.out.counts
            .flatMap { _meta, table, summary -> [ table, summary ] }
            .toSortedList { a, b -> a.name <=> b.name },
        ch_rrna_logs.map { _meta, log -> log }.toSortedList { a, b -> a.name <=> b.name },
        COVERM_GENOME.out.stats.map { _meta, stats -> stats }.ifEmpty([]),
        LINK_CONTIGS.out.links.map { _meta, links -> links },
        ch_all_contig2mag,
        ch_summary,
        ch_gene_clusters.map { _meta, clusters -> clusters }.first(),
        PREPARE_REFERENCE.out.saf.map { _meta, saf -> saf }
    )
    ch_versions = ch_versions.mix(AGGREGATE_COUNTS.out.versions)

    // The tables the database is built from
    ch_database_tables = AGGREGATE_COUNTS.out.libraries
        .mix(
            LINK_CONTIGS.out.links,
            AGGREGATE_COUNTS.out.contig_dna_counts,
            AGGREGATE_COUNTS.out.link_dna_counts,
            AGGREGATE_COUNTS.out.mag_dna_profile,
            AGGREGATE_COUNTS.out.gene_rna_counts
        )
        .map { _meta, table -> table }
        .toSortedList { a, b -> a.name <=> b.name }

    emit:
    contig_links      = LINK_CONTIGS.out.links                 // channel: [ meta, contig_links.tsv ]
    libraries         = AGGREGATE_COUNTS.out.libraries         // channel: [ meta, libraries.tsv ]
    contig_dna_counts = AGGREGATE_COUNTS.out.contig_dna_counts // channel: [ meta, contig_dna_counts.tsv ]
    link_dna_counts   = AGGREGATE_COUNTS.out.link_dna_counts   // channel: [ meta, link_dna_counts.tsv ]
    mag_dna_profile   = AGGREGATE_COUNTS.out.mag_dna_profile   // channel: [ meta, mag_dna_profile.tsv ]
    final_mag_dna     = AGGREGATE_COUNTS.out.final_mag_dna     // channel: [ meta, final_mag_dna.tsv ]
    fine_group_dna    = AGGREGATE_COUNTS.out.fine_group_dna    // channel: [ meta, fine_group_dna.tsv ]
    gene_rna_counts   = AGGREGATE_COUNTS.out.gene_rna_counts   // channel: [ meta, gene_rna_counts.tsv ]
    gene_expression   = AGGREGATE_COUNTS.out.gene_expression   // channel: [ meta, gene_expression.tsv ]
    bams              = ch_bams                                // channel: [ meta, bam, bai ]
    database_tables   = ch_database_tables                     // channel: [ tables of the database ]
    versions          = ch_versions                            // channel: [ versions.yml ]
}
