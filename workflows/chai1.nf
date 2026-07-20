/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Chai1
//
include { RUN_CHAI1 } from '../modules/local/run_chai1'
include { modeChannel } from '../subworkflows/local/utils_nfcore_proteinfold_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHAI1 {

    take:
    ch_samplesheet  // channel: samplesheet read from --input
    ch_versions     // channel: [ path(versions.yml) ]
    ch_chai1_conformers
    ch_chai1_models
    ch_chai1_esm

    main:

    // CHAI1
    RUN_CHAI1(
        ch_samplesheet,
        ch_chai1_conformers,
        ch_chai1_models,
        ch_chai1_esm
    )
    ch_versions = ch_versions.mix(RUN_CHAI1.out.versions)

    RUN_CHAI1
        .out
        .multiqc
        .map { it -> it[1] }
        .toSortedList()
        .map { it ->
            [[ "model": "chai1" ], it.flatten() ]
        }
        .set { ch_multiqc_report }

    modeChannel(RUN_CHAI1.out.top_ranked_pdb, "chai1").set { ch_pdb_final }

    emit:
    versions       = ch_versions
    pdb            = ch_pdb_final
    multiqc_report = ch_multiqc_report
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
