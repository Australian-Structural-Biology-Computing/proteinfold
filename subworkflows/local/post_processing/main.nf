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
include { MODELCIF_ASSEMBLE_BATCH   } from '../../../modules/local/assemble_modelcif/main'
include { MODELCIF_VALIDATE         } from '../../../modules/local/modelcif_validate/main'
include { CIFCHECK                  } from '../../../modules/local/cifcheck/main'

workflow POST_PROCESSING {

    take:
    skip_visualisation
    ch_report_template
    skip_multiqc
    outdir
    ch_multiqc_rep
    ch_multiqc_config
    ch_multiqc_custom_config
    multiqc_logo
    ch_multiqc_methods_description
    ch_software_versions
    write_modelcif
    modelcif_binary
    modelcif_pae_embed
    modelcif_software_details
    ch_packed_batches

    main:
    def ch_multiqc_files = channel.empty()
    ch_modelcif_published = channel.empty()

    // ModelCIF assembly reads the packed metric stores natively; VALIDATE /
    // CIFCHECK then validate and publish each batch's files. Skipped under
    // stubs: the validator parses real mmCIF content and stub placeholders
    // are empty by design.
    if (write_modelcif && !workflow.stubRun) {
        MODELCIF_ASSEMBLE_BATCH(
            ch_packed_batches.combine(ch_software_versions),
            modelcif_software_details
        )
        MODELCIF_VALIDATE(MODELCIF_ASSEMBLE_BATCH.out.modelcif)
        CIFCHECK(
            MODELCIF_VALIDATE.out.modelcif.map { meta, modelcif ->
                [meta, modelcif, file("$projectDir/bin/dicts/mmcif_ma.sdb", checkIfExists: true)]
            }
        )
        ch_modelcif_published = CIFCHECK.out.modelcif
    }

    if (!skip_visualisation){
        // One report task per batch; the report generator reads the packed
        // metric stores natively and renders every sample in the batch.
        GENERATE_REPORT_BATCH(
            ch_packed_batches,
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
        GENERATE_MULTIQC_CONTENTS(
            ch_multiqc_rep,
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
