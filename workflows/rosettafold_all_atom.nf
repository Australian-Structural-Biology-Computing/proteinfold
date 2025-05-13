/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { RUN_ROSETTAFOLD_ALL_ATOM } from '../modules/local/run_rosettafold_all_atom'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ROSETTAFOLD_ALL_ATOM {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    ch_versions             // channel: [ path(versions.yml) ]
    ch_bfd                  // channel: path(bfd)
    ch_uniref30             // channel: path(uniref30)
    ch_pdb100               // channel: path(pdb100)
    ch_rfaa_paper_weights   // channel: path(rfaa_paper_weightsch_dummy_file           // channel: path(NO_file)
    ch_dummy_file           // channel: path(NO_FILE)

    main:
    ch_multiqc_files  = Channel.empty()
    ch_top_ranked_pdb = Channel.empty()
    ch_pdb            = Channel.empty()
    ch_plddt          = Channel.empty()
    ch_msa            = Channel.empty()
    ch_paes           = Channel.empty()
    ch_multiqc_report = Channel.empty()

    RUN_ROSETTAFOLD_ALL_ATOM (
        ch_samplesheet,
        ch_bfd,
        ch_uniref30,
        ch_pdb100,
        ch_rfaa_paper_weights
    )
    ch_versions = ch_versions.mix(RUN_ROSETTAFOLD_ALL_ATOM.out.versions)

    RUN_ROSETTAFOLD_ALL_ATOM
        .out
        .pdb
        .combine(ch_dummy_file)
        .map {
            it[0]["model"] = "rosettafold_all_atom"
            [ it[0]["id"], it[0], it[1], it[2] ]
        }
        .set { ch_top_ranked_pdb }

    RUN_ROSETTAFOLD_ALL_ATOM
        .out
        .plddt
        .map { it[1] }
        .toSortedList()
        .map { [ [ "model": "rosettafold_all_atom" ], it.flatten() ] }
        .set { ch_multiqc_report }

    emit:
    top_ranked_pdb = ch_top_ranked_pdb // channel: [ id, /path/to/*.pdb ]
    pdb            = ch_pdb
    plddt          = ch_plddt
    msa            = ch_msa
    paes           = ch_paes
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
