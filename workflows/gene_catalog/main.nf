#!/usr/bin/env nextflow

include { ASSIGN_GENES                           } from '../../modules/local/assign_genes'
include { BUILD_CATALOG                          } from '../../modules/local/build_catalog'
include { MMSEQS_CLUSTER as MMSEQS_CLUSTER_MAG   } from '../../modules/local/mmseqs/cluster'
include { MMSEQS_CLUSTER as MMSEQS_CLUSTER_OTHER } from '../../modules/local/mmseqs/cluster'
include { MMSEQS_SEARCH                          } from '../../modules/local/mmseqs/search'
include { PARTITION_GENES                        } from '../../modules/local/partition_genes'

//
// One non-redundant gene catalog, built so that a final-MAG gene is only ever represented by
// a final-MAG gene. Representatives cannot be swapped after a plain clustering of every gene:
// MMseqs2 only guarantees that each member meets the thresholds against its own representative,
// so a member put forward in its place could leave the others below them. The catalog is
// seeded with the final MAGs instead:
//
//   1. the final-MAG genes are clustered among themselves
//   2. every other gene is searched against their representatives, and joins its best hit
//   3. the genes left over are clustered among themselves
//
workflow GENE_CATALOG {

    take:
    ch_gff              // channel: [ meta, gff ], one entry per assembly
    ch_fna              // channel: [ meta, fna ], one entry per assembly
    ch_faa              // channel: [ meta, faa ], one entry per assembly
    ch_contig2mag       // channel: [ meta, contig2mag.tsv ], one entry per assembly
    ch_mags_summary     // channel: [ meta, mags_summary.tsv ]
    ch_final_contig2mag // channel: [ meta, final_contig2mag.tsv ]

    main:
    ch_versions = channel.empty()

    //
    // MODULE: split the proteins of each assembly between the final MAGs and the rest
    //
    // `first()` makes the run-wide table a value channel: as the queue channel a process
    // emits, its single item would be paired with the first assembly only.
    PARTITION_GENES(ch_faa, ch_final_contig2mag.map { _meta, final_contig2mag -> final_contig2mag }.first())
    ch_versions = ch_versions.mix(PARTITION_GENES.out.versions.first())

    // Sorted for the same reason as in MAG_DEREPLICATION: so that the run-wide tasks keep
    // their hash, and `-resume` its cache, whatever order the assemblies finish in.
    ch_mag_genes = PARTITION_GENES.out.mag_genes
        .map { _meta, faa -> faa }
        .filter { faa -> faa.size() > 0 }
        .toSortedList { a, b -> a.name <=> b.name }

    ch_other_genes = PARTITION_GENES.out.other_genes
        .map { _meta, faa -> faa }
        .filter { faa -> faa.size() > 0 }
        .toSortedList { a, b -> a.name <=> b.name }

    //
    // MODULE: step 1, cluster the final-MAG genes among themselves (skipped without a final MAG)
    //
    MMSEQS_CLUSTER_MAG(
        ch_mag_genes
            .filter { faas -> faas }
            .map { faas -> [ [ id: 'final_mag_genes' ], faas ] }
    )
    ch_versions = ch_versions.mix(MMSEQS_CLUSTER_MAG.out.versions)

    //
    // MODULE: step 2, search every other gene against the final-MAG representatives
    //
    MMSEQS_SEARCH(
        ch_other_genes
            .filter { faas -> faas }
            .map { faas -> [ [ id: 'other_genes' ], faas ] },
        MMSEQS_CLUSTER_MAG.out.rep_seq.map { _meta, rep_seq -> rep_seq }
    )
    ch_versions = ch_versions.mix(MMSEQS_SEARCH.out.versions)

    //
    // MODULE: hand each matching gene over to its best final-MAG representative
    //
    ASSIGN_GENES(
        ch_other_genes.map { faas -> [ [ id: 'other_genes' ], faas ] },
        MMSEQS_SEARCH.out.hits.map { _meta, hits -> hits }.ifEmpty([])
    )
    ch_versions = ch_versions.mix(ASSIGN_GENES.out.versions)

    //
    // MODULE: step 3, cluster the genes no final-MAG representative took
    //
    MMSEQS_CLUSTER_OTHER(
        ASSIGN_GENES.out.unassigned.filter { _meta, faa -> faa.size() > 0 }
    )
    ch_versions = ch_versions.mix(MMSEQS_CLUSTER_OTHER.out.versions)

    //
    // MODULE: merge the three steps into the catalog
    //
    ch_predictions = ch_gff
        .join(ch_fna)
        .join(ch_faa)
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .map { assemblies ->
            [ [ id: 'gene_catalog' ],
              assemblies.collect { assembly -> assembly[1] },
              assemblies.collect { assembly -> assembly[2] },
              assemblies.collect { assembly -> assembly[3] } ]
        }

    BUILD_CATALOG(
        ch_predictions,
        MMSEQS_CLUSTER_MAG.out.clusters.map { _meta, clusters -> clusters }.ifEmpty([]),
        ASSIGN_GENES.out.assigned.map { _meta, assigned -> assigned },
        MMSEQS_CLUSTER_OTHER.out.clusters.map { _meta, clusters -> clusters }.ifEmpty([]),
        ch_contig2mag
            .map { _meta, contig2mag -> contig2mag }
            .toSortedList { a, b -> a.name <=> b.name },
        ch_mags_summary.map { _meta, summary -> summary }
    )
    ch_versions = ch_versions.mix(BUILD_CATALOG.out.versions)

    emit:
    faa      = BUILD_CATALOG.out.faa      // channel: [ meta, gene_catalog.faa ]
    fna      = BUILD_CATALOG.out.fna      // channel: [ meta, gene_catalog.fna ]
    gff      = BUILD_CATALOG.out.gff      // channel: [ meta, gene_catalog.gff ]
    clusters = BUILD_CATALOG.out.clusters // channel: [ meta, gene_clusters.tsv ]
    stats    = BUILD_CATALOG.out.stats    // channel: [ meta, gene_catalog_stats.tsv ]
    versions = ch_versions                // channel: [ versions.yml ]
}
