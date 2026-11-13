//
// Subworkflow with functionality specific to the dsamoht/clustermg pipeline
//

include { paramsHelp         } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'
include { samplesheetToList  } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    validate_params // boolean: validate parameters against nextflow_schema.json
    input           // string : path to the sample sheet
    reads           // string : path to the sheet of read libraries, or null
    help            // boolean or string: display the help message
    help_full       // boolean: display the full help message
    show_hidden     // boolean: display hidden parameters in the help message
    version         // boolean: display the pipeline version

    main:

    //
    // Banner: same ASCII logo and description as the README, printed on top of the help
    // message and of every run
    //
    log.info(logo())

    //
    // Print the help message built from nextflow_schema.json and exit
    //
    if (help || help_full) {
        log.info(paramsHelp(
            [
                command   : "nextflow run dsamoht/clustermg -profile <docker/singularity/apptainer> --input samplesheet.csv --outdir <OUTDIR>",
                showHidden: show_hidden,
                fullHelp  : help_full,
            ],
            help instanceof String && help != "true" ? help : "",
        ))
        removeEmptyWorkDir()
        exit(0)
    }

    //
    // Print the pipeline version and exit
    //
    if (version) {
        log.info("${workflow.manifest.name} v${workflow.manifest.version}")
        removeEmptyWorkDir()
        exit(0)
    }

    //
    // Validate parameters and print a summary of the ones that differ from the defaults
    //
    if (validate_params) {
        validateParameters()
    }

    // `outputDir` follows `params.outdir` (nextflow.config), and Nextflow falls back to
    // publishing in `./results` when it is null. Parameter validation already requires
    // `--outdir`, but it can be turned off, so the run is stopped here as well: the results
    // are always written to the directory asked for on the command line, never anywhere else.
    if (!params.outdir) {
        error("No output directory: pass --outdir. Without it the results would be published to './results' instead.")
    }

    // The finer MAG groups are cut from dRep's tree below its own threshold, never above
    if (params.mag_ani_fine < params.mag_ani) {
        error("--mag_ani_fine (${params.mag_ani_fine}) must be at least --mag_ani (${params.mag_ani}): the finer groups are cut inside dRep's clusters.")
    }

    log.info(paramsSummaryLog(workflow))

    //
    // Read the sample sheet. Per-row validation lives in assets/schema_input.json;
    // nf-schema returns [ meta, assembly, bins, checkm2, samples ] as strings, with an empty
    // list for a column left blank. The paths are resolved here rather than by nf-schema, so
    // that a relative path can also be found next to the sample sheet.
    //
    def samplesheet_dir = file(input).parent
    def sheet = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, assembly, bins, checkm2, sample_ids ->
            [
                meta,
                resolveInputPath(assembly, 'assembly', meta.id, samplesheet_dir, false),
                bins    ? resolveInputPath(bins, 'bins', meta.id, samplesheet_dir, true)     : [],
                checkm2 ? resolveInputPath(checkm2, 'checkm2', meta.id, samplesheet_dir, false) : [],
                sample_ids ? sample_ids.tokenize(';') : []
            ]
        }
    def samples = sheet.collect { meta, assembly, bins, checkm2, _sample_ids -> [ meta, assembly, bins, checkm2 ] }

    //
    // Checks that span several rows of the sample sheet, which a per-row JSON schema
    // cannot express: the schema only sees paths as written, two spellings of the same
    // file only meet once resolved. The same assembly twice would put every one of its
    // genes and MAGs in the run twice.
    //
    samples
        .groupBy { _meta, assembly, _bins, _checkm2 -> assembly.toUriString() }
        .findAll { _path, rows -> rows.size() > 1 }
        .each { path, rows ->
            error("assembly '${path}' is given for several rows: ${rows.collect { row -> row[0].id }.join(', ')}. Each assembly must appear once.")
        }

    ch_samplesheet = channel.fromList(samples)

    //
    // Every assembly, the samples it was built from and whether that makes it a single-sample
    // assembly or a co-assembly: provenance for the database, it changes no result
    //
    def assembly_rows = sheet.collect { meta, assembly, _bins, _checkm2, sample_ids ->
        def kind = sample_ids.size() > 1 ? 'co' : (sample_ids ? 'single' : '')
        [ meta.id, kind, sample_ids.join(';'), assembly.toUriString() ].join('\t')
    }
    ch_assemblies = channel.of(([ 'assembly\tkind\tsamples\tsource' ] + assembly_rows).join('\n') + '\n')
        .collectFile(name: 'assemblies.tsv')

    //
    // Read the sheet of read libraries, validated by assets/schema_reads.json. nf-schema
    // returns [ meta, fastq_1, fastq_2 ], the library id, sample, molecule, platform and
    // strandedness in meta.
    //
    def libraries = []
    if (reads) {
        def reads_dir = file(reads).parent
        libraries = samplesheetToList(reads, "${projectDir}/assets/schema_reads.json")
            .collect { meta, fastq_1, fastq_2 ->
                def files = [ resolveInputPath(fastq_1, 'fastq_1', meta.id, reads_dir, false) ]
                if (fastq_2) {
                    files << resolveInputPath(fastq_2, 'fastq_2', meta.id, reads_dir, false)
                }
                [ meta + [ strandedness: meta.strandedness ?: null, single_end: files.size() == 1 ], files ]
            }
    }

    //
    // What the schema cannot express on its own. RNA is counted on genes from short reads,
    // with its strand; DNA is counted on contigs, strand aside.
    //
    libraries.each { meta, files ->
        if (meta.molecule == 'rna' && meta.platform != 'illumina') {
            error("library '${meta.id}': RNA is counted on genes from Illumina reads; platform '${meta.platform}' is not supported for RNA.")
        }
        if (meta.molecule == 'rna' && !meta.strandedness) {
            error("library '${meta.id}': an RNA library needs its 'strandedness': unstranded, forward or reverse.")
        }
        if (meta.molecule == 'dna' && meta.strandedness) {
            error("library '${meta.id}': 'strandedness' only applies to RNA libraries.")
        }
        if (files.size() == 2 && meta.platform != 'illumina') {
            error("library '${meta.id}': only Illumina libraries have mates ('fastq_2').")
        }
    }

    // A FastQ file given twice would have its reads mapped, and counted, twice
    libraries
        .collectMany { meta, files -> files.collect { f -> [ f.toUriString(), meta.id ] } }
        .groupBy { uri, _id -> uri }
        .findAll { _uri, entries -> entries.size() > 1 }
        .each { uri, entries ->
            error("FastQ '${uri}' is given more than once (${entries.collect { entry -> entry[1] }.join(', ')}): its reads would be counted twice.")
        }

    // The samples an assembly was built from are expected among the libraries, but a run
    // may well leave some out
    if (libraries) {
        def library_samples = libraries.collect { meta, _files -> meta.sample }.toSet()
        def unmapped = sheet.collectMany { row -> row[4] }.unique().findAll { sample_id -> !(sample_id in library_samples) }
        if (unmapped) {
            log.warn("No read library for sample(s) ${unmapped.join(', ')}, listed in the 'samples' of the assemblies.")
        }
    }

    ch_reads = channel.fromList(libraries)

    //
    // The rRNA databases SortMeRNA filters the RNA libraries against, read from their manifest
    // only when there is RNA to filter
    //
    def rna = libraries.any { meta, _files -> meta.molecule == 'rna' }
    def rrna_fastas = rna && !params.skip_rrna_removal ? riboDatabases(params.ribo_database_manifest) : []
    ch_rrna_databases = channel.value(rrna_fastas)

    emit:
    samplesheet    = ch_samplesheet
    assemblies     = ch_assemblies     // channel: assemblies.tsv, the kind and samples of every assembly
    reads          = ch_reads          // channel: [ meta, [ fastq ] ], one entry per read library
    rrna_databases = ch_rrna_databases // channel: [ rRNA FASTA ], empty without RNA to filter
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// ASCII logo, pipeline description and repository details, taken from the manifest so
// that the banner cannot drift from the version being run.
//
def logo() {
    return """

         __              __
  _____ / /__  __ _____ / /_ ___   _____ ____ ___   ____ _
 / ___// // / / // ___// __// _ \\ / ___// __ `__ \\ / __ `/
/ /__ / // /_/ /(__  )/ /_ /  __// /   / / / / / // /_/ /
\\___//_/ \\__,_//____/ \\__/ \\___//_/   /_/ /_/ /_/ \\__, /
                                                 /____/

${workflow.manifest.description}

     Github: ${workflow.manifest.homePage}
     Version: v${workflow.manifest.version}

"""
}

//
// Nextflow creates the work directory when the session starts, before this script runs, so
// a run that only prints the help message or the version leaves an empty `work/` behind in
// whatever directory it was called from. Remove it again, but only if it is still empty,
// i.e. if this run is the one that created it.
//
def removeEmptyWorkDir() {
    if (workflow.workDir.exists() && workflow.workDir.list().size() == 0) {
        workflow.workDir.deleteDir()
    }
}

//
// Resolve a path of the sample sheet. A relative path is looked up in the launch directory
// first, as nf-schema would, then next to the sample sheet, so that a sample sheet can be
// kept with its data and used from anywhere. Remote paths (`s3://`, `https://`, ...) are
// taken as they are.
//
def resolveInputPath(String value, String column, String name, samplesheet_dir, Boolean is_dir) {
    def path = file(value)
    def relative = !(value ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/) && !value.startsWith('/')
    if (!path.exists() && relative && samplesheet_dir.resolve(value).exists()) {
        path = samplesheet_dir.resolve(value)
    }
    if (!path.exists()) {
        def tried = relative ? " (looked in the launch directory and in ${samplesheet_dir})" : ''
        error("row '${name}': ${column} '${value}' does not exist${tried}")
    }
    if (is_dir && !path.isDirectory()) {
        error("row '${name}': ${column} '${value}' is not a directory")
    }
    if (!is_dir && path.isDirectory()) {
        error("row '${name}': ${column} '${value}' is a directory, not a file")
    }
    return path
}

//
// Read the manifest of rRNA databases: one FASTA per line, a path or a URL. A relative path is
// looked up next to the manifest; blank lines and lines starting with '#' are skipped.
//
def riboDatabases(String manifest) {
    def manifest_file = file(manifest)
    if (!manifest_file.exists()) {
        error("rRNA database manifest '${manifest}' does not exist")
    }
    def databases = manifest_file.readLines()
        .collect { line -> line.trim() }
        .findAll { line -> line && !line.startsWith('#') }
        .collect { line ->
            def remote = line ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/
            def path = remote || line.startsWith('/') ? file(line) : manifest_file.parent.resolve(line)
            if (!remote && !path.exists()) {
                error("rRNA database '${line}' of manifest '${manifest}' does not exist")
            }
            path
        }
    if (!databases) {
        error("rRNA database manifest '${manifest}' lists no FASTA file")
    }
    return databases
}
