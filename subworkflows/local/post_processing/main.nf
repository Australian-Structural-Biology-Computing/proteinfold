//
// Post processing analysis for the predicted structures
//

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../utils_nfcore_proteinfold_pipeline'

include { GENERATE_REPORT_BATCH    } from '../../../modules/local/generate_report'
include { GENERATE_MULTIQC_CONTENTS } from '../../../modules/local/generate_multiqc_contents'
include { MULTIQC                   } from '../../../modules/nf-core/multiqc/main'
include { MODELCIF_ASSEMBLE         } from '../../../modules/local/assemble_modelcif/main'
include { MODELCIF_VALIDATE         } from '../../../modules/local/modelcif_validate/main'
include { CIFCHECK                  } from '../../../modules/local/cifcheck/main'

workflow POST_PROCESSING {

    take:
    skip_visualisation
    ch_report_input
    ch_report_template
    skip_multiqc
    outdir
    ch_multiqc_rep
    ch_multiqc_config
    ch_multiqc_custom_config
    multiqc_logo
    ch_multiqc_methods_description
    ch_software_versions
    ch_modelcif
    write_modelcif
    modelcif_binary
    modelcif_pae_embed
    analysis_batch_size // int: number of samples analysed together (-1 = one batch)

    main:
    if (analysis_batch_size == 0 || analysis_batch_size < -1) {
        error("Analysis batch size must be -1 or a positive integer, got ${analysis_batch_size}")
    }
    def ch_multiqc_files = channel.empty()
    ch_modelcif_published = channel.empty()

    // ModelCIF assembly consumes real ranked-structure metrics; stub runs
    // only produce empty placeholder files, so skip assembly entirely there.
    if (write_modelcif && !workflow.stubRun) {
        MODELCIF_ASSEMBLE(ch_modelcif)
        MODELCIF_VALIDATE(MODELCIF_ASSEMBLE.out.modelcif)
        CIFCHECK(
            MODELCIF_VALIDATE.out.modelcif.map { meta, modelcif ->
                [meta, modelcif, file("$projectDir/bin/dicts/mmcif_ma.sdb", checkIfExists: true)]
            }
        )
        ch_modelcif_published = CIFCHECK.out.modelcif
    }

    if (!skip_visualisation){
        // Group the per-sample report inputs by public mode (the generator's
        // --type/--tag are constant within a mode) and collate each mode's
        // samples into --analysis_batch_size-sized report tasks. The key is
        // the full mode name: distinct modes can share a mode_short token,
        // and their per-sample files would collide inside one task.
        def ch_report_batches = ch_report_input
            .map { meta, structure, msa, pae, iptm, ipsae, chainwise, chainwise_ipsae ->
                tuple(meta.model, meta, structure, msa, pae, iptm, ipsae, chainwise, chainwise_ipsae)
            }
            .groupTuple()
            .flatMap { mode, metas, structures, msas, paes, iptms, ipsaes, chainwises, chainwise_ipsaes ->
                def rows = new ArrayList([metas, structures, msas, paes, iptms, ipsaes, chainwises, chainwise_ipsaes])
                    .transpose()
                    .sort { a, b -> (a[0].id as String) <=> (b[0].id as String) }
                def effective = analysis_batch_size == -1 ? Math.max(1, rows.size()) : analysis_batch_size
                rows.collate(effective).collect { batch ->
                    tuple(
                        new ArrayList(batch*.getAt(0)),
                        // The ranked-structure slot is a per-sample list of
                        // lists; the path input needs one flat file list.
                        new ArrayList(batch*.getAt(1).flatten()),
                        new ArrayList(batch*.getAt(2)),
                        new ArrayList(batch*.getAt(3)),
                        new ArrayList(batch*.getAt(4)),
                        new ArrayList(batch*.getAt(5)),
                        new ArrayList(batch*.getAt(6)),
                        new ArrayList(batch*.getAt(7))
                    )
                }
            }
        GENERATE_REPORT_BATCH(
            ch_report_batches,
            ch_report_template
        )
    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = channel.empty()

    if (!skip_multiqc) {
        ch_summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary    = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files       = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_methods_description))
        ch_multiqc_files       = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        ch_multiqc_files       = ch_multiqc_files.mix(ch_software_versions)

        def mqc_generator_script = file("${projectDir}/modules/local/generate_multiqc_contents/generate_multiqc_contents.py", checkIfExists: true)
        // Chunk each mode's metrics so one task never has to hold every
        // sample's per-residue series; unsplit modes keep a single task with
        // the unchunked section ids.
        def ch_multiqc_chunks = ch_multiqc_rep
            .flatMap { meta, files ->
                def list = files instanceof List ? files : [files]
                def effective = analysis_batch_size == -1 ? Math.max(1, list.size()) : analysis_batch_size
                def chunks = new ArrayList(list).sort { it.name }.collate(effective)
                if (chunks.size() <= 1) {
                    // flatMap flattens one level: an unsplit mode must come
                    // back as a single tuple, not as its two elements.
                    return [[meta, list]]
                }
                chunks.withIndex().collect { chunk, idx -> [meta + [analysis_chunk: idx], new ArrayList(chunk)] }
            }
        GENERATE_MULTIQC_CONTENTS(
            ch_multiqc_chunks,
            mqc_generator_script
        )

        // Multi-mode runs emit one contents task per mode; collect them so
        // MULTIQC runs once over every mode's custom content.
        MULTIQC (
            GENERATE_MULTIQC_CONTENTS.out.mqc_json
                .map { meta, files -> [ [id: 'multiqc', model: 'proteinfold'], files instanceof List ? files : [files] ] }
                .groupTuple()
                .map { meta, files -> [ meta, files.flatten() ] }
                .combine(ch_multiqc_files.collect().map { files -> [files] })
                .combine(ch_multiqc_config.collect().ifEmpty([]).map { config -> [config] })
                .combine(ch_multiqc_custom_config.collect().ifEmpty([]).map { config -> [config] })
                .map { meta, mqc_json, extra_files, config_files, custom_config_files ->
                    // A single-file glob output arrives as a bare java.nio.file.Path, which is Iterable: `path + list` would splice its name components into the list.
                    def mqc_json_files = mqc_json instanceof List ? mqc_json : [mqc_json]
                    [
                        meta,
                        mqc_json_files + extra_files,
                        (config_files ?: []) + (custom_config_files ?: []),
                        multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                        [],
                        []
                    ]
                }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> report }
    }

    emit:
    multiqc_report = ch_multiqc_report
    modelcif = ch_modelcif_published
}
