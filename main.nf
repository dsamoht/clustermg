#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dsamoht/clustermg
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/dsamoht/clustermg
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_DATABASE          } from './modules/local/build_database'
include { GENE_CATALOG            } from './workflows/gene_catalog'
include { MAG_DEREPLICATION       } from './workflows/mag_dereplication'
include { MERGE_PYRODIGAL         } from './modules/local/pyrodigal/merge'
include { PIPELINE_INITIALISATION } from './subworkflows/local/pipeline_initialisation'
include { PREPARE_ASSEMBLY        } from './modules/local/prepare_assembly'
include { PYRODIGAL               } from './modules/local/pyrodigal/pyrodigal'
include { READ_MAPPING            } from './workflows/read_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOOLEAN PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Declared with their type so that Nextflow converts `--version false` on the command line
// to false: undeclared, the value stays the string "false", which Groovy reads as true.
// Defaults are set in nextflow.config. `--help` is left out, it also takes a string.
params {
    monochrome_logs : Boolean
    help_full       : Boolean
    show_hidden     : Boolean
    version         : Boolean
    validate_params   : Boolean
    save_bams         : Boolean
    skip_rrna_removal : Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLUSTERMG {

    take:
    ch_samplesheet // channel: [ meta, assembly, bins, checkm2 ], from the sample sheet
    ch_assemblies  // channel: assemblies.tsv, the kind and samples of every assembly
    ch_reads       // channel: [ meta, [ fastq ] ], one entry per read library, none without --reads
    ch_rrna        // channel: [ rRNA FASTA ], value, the databases the RNA reads are filtered against

    main:
    ch_versions = channel.empty()

    //
    // MODULE: rename the contigs and MAGs of each assembly, check the bins against it
    //
    PREPARE_ASSEMBLY(ch_samplesheet)
    ch_versions = ch_versions.mix(PREPARE_ASSEMBLY.out.versions.first())

    //
    // MODULE: predict genes on 100 MB chunks of each assembly
    //
    ch_fasta_chunks = PREPARE_ASSEMBLY.out.contigs
        .map { meta, contigs -> [ meta.id, contigs ] }
        .splitFasta(size: 100.MB, file: true)
        .map { name, chunk -> [ name, chunk.baseName.tokenize('.').last(), chunk ] }

    // Number of chunks per assembly, so MERGE_PYRODIGAL can start on an assembly as soon
    // as that assembly's own chunks are annotated
    ch_chunks_per_assembly = ch_fasta_chunks
        .map { name, _chunk_idx, _chunk -> [ name, 1 ] }
        .groupTuple()
        .map { name, chunks -> [ name, chunks.size() ] }

    ch_gene_chunks = PYRODIGAL(
        ch_fasta_chunks.map { name, chunk_idx, chunk -> [ [ id: name ], chunk_idx, chunk ] }
    )
    ch_versions = ch_versions.mix(PYRODIGAL.out.versions.first())

    // Pyrodigal emits the gff/fna/faa of a chunk as a single tuple, so they stay paired,
    // and the chunk index travels with them instead of being parsed back from the file
    // names at merge time.
    ch_genes_to_merge = ch_gene_chunks.genes
        .map { meta, chunk_idx, gff, fna, faa -> [ meta.id, chunk_idx, gff, fna, faa ] }
        .combine(ch_chunks_per_assembly, by: 0)
        .map { name, chunk_idx, gff, fna, faa, n_chunks ->
            [ groupKey(name, n_chunks), chunk_idx as Integer, gff, fna, faa ]
        }
        .groupTuple()
        .map { name, chunk_idxs, gffs, fnas, faas ->
            def ordered = [ chunk_idxs, gffs, fnas, faas ].transpose().sort { a, b -> a[0] <=> b[0] }
            [ [ id: name.toString() ],
              ordered.collect { chunk -> chunk[1] },
              ordered.collect { chunk -> chunk[2] },
              ordered.collect { chunk -> chunk[3] } ]
        }

    ch_genes = MERGE_PYRODIGAL(ch_genes_to_merge)
    ch_versions = ch_versions.mix(MERGE_PYRODIGAL.out.versions.first())

    //
    // SUBWORKFLOW: dereplicate the MAGs of every assembly into one final set
    //
    MAG_DEREPLICATION(
        PREPARE_ASSEMBLY.out.mags,
        PREPARE_ASSEMBLY.out.mag_info,
        PREPARE_ASSEMBLY.out.contig2mag
    )
    ch_versions = ch_versions.mix(MAG_DEREPLICATION.out.versions)

    //
    // SUBWORKFLOW: build the gene catalog around the final MAGs
    //
    GENE_CATALOG(
        ch_genes.gff,
        ch_genes.fna,
        ch_genes.faa,
        PREPARE_ASSEMBLY.out.contig2mag,
        MAG_DEREPLICATION.out.summary,
        MAG_DEREPLICATION.out.final_contig2mag
    )
    ch_versions = ch_versions.mix(GENE_CATALOG.out.versions)

    //
    // SUBWORKFLOW: map every read library once onto all the assemblies, and count it
    //
    READ_MAPPING(
        ch_reads,
        PREPARE_ASSEMBLY.out.contigs,
        ch_genes.gff,
        PREPARE_ASSEMBLY.out.contig2mag,
        PREPARE_ASSEMBLY.out.mags,
        MAG_DEREPLICATION.out.summary,
        GENE_CATALOG.out.clusters,
        ch_rrna
    )
    ch_versions = ch_versions.mix(READ_MAPPING.out.versions)

    //
    // MODULE: load every result into one DuckDB database
    //
    BUILD_DATABASE(
        ch_assemblies.map { assemblies -> [ [ id: 'clustermg' ], assemblies ] },
        PREPARE_ASSEMBLY.out.contig_names.map { _meta, names -> names }.toSortedList { a, b -> a.name <=> b.name },
        PREPARE_ASSEMBLY.out.contig2mag.map { _meta, contig2mag -> contig2mag }.toSortedList { a, b -> a.name <=> b.name },
        MAG_DEREPLICATION.out.summary.map { _meta, summary -> summary },
        ch_genes.gff.map { _meta, gff -> gff }.toSortedList { a, b -> a.name <=> b.name },
        GENE_CATALOG.out.clusters.map { _meta, clusters -> clusters },
        GENE_CATALOG.out.faa.map { _meta, faa -> faa },
        GENE_CATALOG.out.fna.map { _meta, fna -> fna },
        READ_MAPPING.out.database_tables
    )
    ch_versions = ch_versions.mix(BUILD_DATABASE.out.versions)

    //
    // Collate software versions
    //
    // Collect the contents rather than the paths: `sort: true` orders the entries by their
    // value, and a channel of paths sorts by work directory, which is a fresh random hash on
    // every run. Reading each file first makes both the deduplication and the order depend on
    // the versions themselves, so the published file is the same from one run to the next.
    ch_collated_versions = ch_versions
        .map { versions -> versions.text }
        .unique()
        .collectFile(name: 'clustermg_software_versions.yml', sort: true)

    emit:
    versions         = ch_collated_versions
    contig_names     = PREPARE_ASSEMBLY.out.contig_names
    predicted_gff    = ch_genes.gff
    predicted_fna    = ch_genes.fna
    predicted_faa    = ch_genes.faa
    final_mags       = MAG_DEREPLICATION.out.final_mags
    mags_summary     = MAG_DEREPLICATION.out.summary
    final_contig2mag = MAG_DEREPLICATION.out.final_contig2mag
    drep_tables      = MAG_DEREPLICATION.out.drep_tables
    catalog_faa      = GENE_CATALOG.out.faa
    catalog_fna      = GENE_CATALOG.out.fna
    catalog_gff      = GENE_CATALOG.out.gff
    gene_clusters    = GENE_CATALOG.out.clusters
    catalog_stats    = GENE_CATALOG.out.stats
    contig2mag       = PREPARE_ASSEMBLY.out.contig2mag
    contig_links     = READ_MAPPING.out.contig_links
    libraries        = READ_MAPPING.out.libraries
    contig_counts    = READ_MAPPING.out.contig_dna_counts.mix(READ_MAPPING.out.link_dna_counts)
    mag_counts       = READ_MAPPING.out.mag_dna_profile.mix(READ_MAPPING.out.final_mag_dna, READ_MAPPING.out.fine_group_dna)
    gene_counts      = READ_MAPPING.out.gene_rna_counts.mix(READ_MAPPING.out.gene_expression)
    bams             = params.save_bams ? READ_MAPPING.out.bams : channel.empty()
    database         = BUILD_DATABASE.out.database.mix(BUILD_DATABASE.out.parquet, BUILD_DATABASE.out.checks)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN THE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: validate the parameters and read the sample sheet
    //
    PIPELINE_INITIALISATION(
        params.validate_params,
        params.input,
        params.reads,
        params.help,
        params.help_full,
        params.show_hidden,
        params.version
    )

    //
    // WORKFLOW: run the main analysis
    //
    CLUSTERMG(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.assemblies,
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.rrna_databases
    )

    publish:
    versions         = CLUSTERMG.out.versions
    contig_names     = CLUSTERMG.out.contig_names
    predicted_gff    = CLUSTERMG.out.predicted_gff
    predicted_fna    = CLUSTERMG.out.predicted_fna
    predicted_faa    = CLUSTERMG.out.predicted_faa
    final_mags       = CLUSTERMG.out.final_mags
    mags_summary     = CLUSTERMG.out.mags_summary
    final_contig2mag = CLUSTERMG.out.final_contig2mag
    drep_tables      = CLUSTERMG.out.drep_tables
    catalog_faa      = CLUSTERMG.out.catalog_faa
    catalog_fna      = CLUSTERMG.out.catalog_fna
    catalog_gff      = CLUSTERMG.out.catalog_gff
    gene_clusters    = CLUSTERMG.out.gene_clusters
    catalog_stats    = CLUSTERMG.out.catalog_stats
    contig2mag       = CLUSTERMG.out.contig2mag
    contig_links     = CLUSTERMG.out.contig_links
    libraries        = CLUSTERMG.out.libraries
    contig_counts    = CLUSTERMG.out.contig_counts
    mag_counts       = CLUSTERMG.out.mag_counts
    gene_counts      = CLUSTERMG.out.gene_counts
    bams             = CLUSTERMG.out.bams
    database         = CLUSTERMG.out.database
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PUBLISH TARGETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    contig_names {
        path "assemblies"
    }
    contig2mag {
        path "assemblies"
    }
    predicted_gff {
        path "genes/predicted"
    }
    predicted_fna {
        path "genes/predicted"
    }
    predicted_faa {
        path "genes/predicted"
    }
    final_mags {
        path "mags/final"
    }
    mags_summary {
        path "mags"
    }
    final_contig2mag {
        path "mags"
    }
    drep_tables {
        path "mags/drep"
    }
    catalog_faa {
        path "genes"
    }
    catalog_fna {
        path "genes"
    }
    catalog_gff {
        path "genes"
    }
    gene_clusters {
        path "genes"
    }
    catalog_stats {
        path "genes"
    }
    contig_links {
        path "coverage"
    }
    libraries {
        path "coverage"
    }
    contig_counts {
        path "coverage"
    }
    mag_counts {
        path "mags"
    }
    gene_counts {
        path "genes"
    }
    bams {
        path "mapping"
    }
    database {
        path "database"
    }
    versions {
        path "pipeline_info"
    }
}
