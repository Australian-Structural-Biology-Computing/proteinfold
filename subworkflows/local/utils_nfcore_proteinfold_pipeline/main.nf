//
// Subworkflow with functionality specific to the nf-core/proteinfold pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'

def sanitise(id) {
    (id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def colors = logColours(monochrome_logs)

    def before_text = ""
    def after_text = ""
    before_text = """
-${colors.dim}----------------------------------------------------${colors.reset}-
                                        ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                        ${colors.green}`._,._,\'${colors.reset}
${colors.purple}  nf-core/proteinfold ${workflow.manifest.version}${colors.reset}
-${colors.dim}----------------------------------------------------${colors.reset}-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/proteinfold/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    validateInputParameters()

    //
    // Create channel from input file provided through input
    //
    def sample_entries = samplesheetToList(input, "assets/schema_input.json")
    def sample_ids = sample_entries.collect { item ->
        def meta = item[0]
        meta.id ?: meta.sequence
    }
    if (sample_ids.size() != sample_ids.toSet().size()) {
        error("Duplicate sample IDs are not allowed: ${sample_ids}")
    }
    ch_samplesheet = channel.fromList(sample_entries)

    ch_samplesheet = ch_samplesheet
        .map { item ->
            def meta = item[0]
            def fasta = item[1]
            def msa = item.size() > 2 ? item[2] : null
            // This mapping supports legacy samplesheets that use 'sequence' as metadata.
            // If meta.id is missing or empty, meta.sequence is used as the identifier.
            def identifier = meta.id ? meta.id : meta.sequence
            def normalized_meta = [id: identifier]
            if (msa) {
                if (params.split_fasta) {
                    error("The 'msa' samplesheet column cannot be used with --split_fasta.")
                }
                normalized_meta.msa = msa
            }
            return [normalized_meta, fasta]
        }

    if (params.split_fasta) {
        ch_samplesheet = ch_samplesheet.map { meta, fasta ->
            validateFasta(fasta)
            [meta, fasta]
        }

        // Split the fasta file into individual files for each sequence
        ch_samplesheet
            .map { meta, fasta -> [meta.id, fasta] }
            .splitFasta( record: [header: true, sequence: true], elem: 1 )
            .collectFile { sample_id, item ->
                [ "${cleanHeader(sample_id.toString())}_${cleanHeader(item["header"])}.fa", ">" + cleanHeader(item["header"]) + '\n' +item["sequence"] ]
            }
            .map {
                file -> [[id: file.baseName], file]
            }
            .set { ch_samplesheet }
    }

    emit:
    samplesheet  = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete = {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError = {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    def requestedModes = params.mode.toLowerCase().split(",").collect { mode -> mode.trim() }

    if (requestedModes.contains("alphafold3")) {
        alphafold3Warn(log)
    }

    if (params.random_seed != null) {
        randomSeedModeWarn(log, requestedModes)
    }
}

def modeChannel(ch, mode) {
    return ch.map { meta, value ->
        def meta_clone = meta.clone()
        meta_clone.model = mode
        [ meta_clone, value ]
    }
}

//
// Collect the per-model metric TSVs that the MultiQC custom content consumes.
// `metric_channels` is a List of [name, channel] pairs. Each channel yields a
// tuple whose [1] and [2] fields are the meta map and the metric file; pinning
// that position is what makes this safe for emits with extra trailing fields
// (boltz emits 5, alphafold2 4). The first non-empty channel seeds the fold and
// the rest are mixed in one at a time, because mix() is a channel operator -- it
// does not exist on ArrayList.
//
// Summarise a batch of sample metas into a bounded task tag: small batches
// keep the full id list, large ones show the first three ids plus the
// remaining count so a 1000-sample batch does not become a 10 kB tag.
def batchTag(batch_metas) {
    def ids = (batch_metas instanceof Collection ? batch_metas.collect() : [batch_metas])
        .collect { (it.id as String) }
    if (ids.size() <= 3) {
        return "batch_${ids.join('_')}"
    }
    return "batch_${ids[0]}_${ids[1]}_${ids[2]}+${ids.size() - 3}more"
}

def batchMergeInputs(ch_merge_input, analysis_batch) {
    // rows: [meta, searched MSA JSON, original input JSON]. One task per
    // analysis batch merges every row; samples that share an input file
    // stage it once and address it through the manifest.
    if (analysis_batch == 0 || analysis_batch < -1) {
        error("Analysis batch size must be -1 or a positive integer, got ${analysis_batch}")
    }
    return ch_merge_input
        .collect(flat: false)
        .filter { rows -> !rows.isEmpty() }
        .flatMap { rows ->
            def sorted = new ArrayList(rows).sort { a, b -> (a[0].id as String) <=> (b[0].id as String) }
            def effective = analysis_batch == -1 ? Math.max(1, sorted.size()) : analysis_batch
            sorted.collate(effective).collect { batch ->
                def originals = batch.collect { it[2] }.unique { it.toString() }
                if (originals.collect { it.name }.toSet().size() != originals.size()) {
                    error("Samples in one merge batch share input file names; rename them or run with --analysis_batch_size 1")
                }
                def manifest = batch.collect { meta, msa_json, original_json ->
                    "${sanitise(meta.id)}\t${msa_json.name}\t${original_json.name}"
                }.join('\n')
                [new ArrayList(batch.collect { it[0] }), new ArrayList(batch.collect { it[1] }), originals, manifest]
            }
        }
}

def collectMultiqcMetrics(metric_channels) {
    // The metric tuples carry their own meta, whose model field names the
    // public mode; grouping collapses every sample of a mode into one
    // tuple keyed by that mode, which is what MULTIQC consumes.
    def acc = null
    metric_channels.each { entry ->
        // Index form: some metric emits carry trailing fields (e.g. the MSA
        // coverage plot), so positional destructuring would fail.
        def ch = entry[1].map { it -> [it[0].model, it[1]] }
        acc = acc == null ? ch : acc.mix(ch)
    }
    if (acc == null) {
        return channel.empty()
    }
    return acc
        .unique { entry -> entry[1] }
        .groupTuple()
        .map { model, paths ->
            [ [id: model, model: model], paths.flatten() ]
        }
}

// Join predictor outputs into the ModelCIF assembly input. Every backend
// emits the PAE/PTM/IPTM metrics; the MSA is the one input that can be
// absent (the ESMFold2 language-model variants), which the inference task
// signals with a zero-size MSA TSV that lands in missing_metrics.
def modelcifInput(structs, msa, plddt, pae, ptm, iptm, versions, softwareDetails, msaTool) {
    structs
        .join(msa)
        .join(plddt)
        .join(pae)
        .join(ptm)
        .join(iptm)
        .combine(versions)
        .map { values ->
            def modelcifMeta = values[0].clone()
            modelcifMeta.msa_tool = msaTool
            modelcifMeta.missing_metrics = [java.nio.file.Files.size(values[2]) == 0 ? 'msa' : null].findAll()
            [modelcifMeta] + values[1..-1] + [softwareDetails]
        }
}

def countMolecularEntitiesInFasta(fasta) {
    return fasta.text
        .readLines()
        .count { line -> line.trim().startsWith('>') }
}

def resolveModelPresetByFastaEntities(fasta, monomerPreset, multimerPreset = 'multimer') {
    return countMolecularEntitiesInFasta(fasta) > 1 ? multimerPreset : monomerPreset
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()
    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def cleanHeader(header) {
    return header
        .replaceAll(" ", "_")
        .replaceAll("/","_")
        .replaceAll(",", "")
        .replaceAll(";","")
        .replaceAll("\\|","_")
}

def validateFasta(fasta) {
    def lines = fasta.text.readLines()
    // extract headers
    def headers = lines.findAll { it -> it.startsWith('>') }
    // if headers are not unique, throw an error
    if (headers.size() != headers.unique().size()) {
        throw new Exception("Invalid FASTA file. The headers are not unique.")
    }
    // check headers that are malformed
    headers.each { header ->
        if (header =~ /[ \t;,\/]/) {
            // warn user that the header contains special characters
            log.warn "The header ${header} contains special characters. They have been automatically removed."
        }
    }
}

//
// Print a warning when using Alphafold3
//
def alphafold3Warn(log) {
    log.warn "=============================================================================\n" +
        "  You are using AlphaFold3 mode.\n" +
        "  Be aware that the predicted structures can not be used for commercial purposes.\n" +
        "  More information here: \"https://github.com/google-deepmind/alphafold3/blob/main/README.md#alphafold-3-source-code-and-model-parameters.\"\n" +
        "==================================================================================="
}

def randomSeedModeWarn(log, requestedModes) {
    // The seed is forwarded to colabfold_batch (--random-seed) and is consumed
    // by the AlphaFold-family prediction backends; the other backends ignore it.
    def supportedModes = [
        'alphafold3',
        'alphafold2'
    ] as Set
    def unsupportedModes = requestedModes.findAll { mode -> !supportedModes.contains(mode) }.unique().sort()

    if (!unsupportedModes.isEmpty()) {
        log.warn "===================================================================================\n" +
            "  The global --random_seed parameter only applies to the AlphaFold-family prediction backends\n" +
            "  (alphafold3 and the alphafold2 modes).\n" +
            "  It will be ignored for the following selected mode(s): ${unsupportedModes.join(', ')}.\n" +
            "==================================================================================="
    }
}
