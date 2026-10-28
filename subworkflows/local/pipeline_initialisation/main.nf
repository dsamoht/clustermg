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

    log.info(paramsSummaryLog(workflow))

    //
    // Read the sample sheet. Per-row validation lives in assets/schema_input.json;
    // nf-schema returns [ meta, assembly, bins, checkm2 ] as strings, with an empty list
    // for a column left blank. The paths are resolved here rather than by nf-schema, so
    // that a relative path can also be found next to the sample sheet.
    //
    def samplesheet_dir = file(input).parent
    def samples = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, assembly, bins, checkm2 ->
            [
                meta,
                resolveInputPath(assembly, 'assembly', meta.id, samplesheet_dir, false),
                bins    ? resolveInputPath(bins, 'bins', meta.id, samplesheet_dir, true)     : [],
                checkm2 ? resolveInputPath(checkm2, 'checkm2', meta.id, samplesheet_dir, false) : []
            ]
        }

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

    emit:
    samplesheet = ch_samplesheet
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
