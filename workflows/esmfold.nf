/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { RUN_ESMFOLD               } from '../modules/local/run_esmfold'
include { MULTIFASTA_TO_SINGLEFASTA } from '../modules/local/multifasta_to_singlefasta'

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

workflow ESMFOLD {

    take:
    ch_samplesheet    // channel: samplesheet read in from --input
    ch_versions       // channel: [ path(versions.yml) ]
    ch_esmfold_params // directory: /path/to/esmfold/params/
    ch_num_recycles   // int: Number of recycles for esmfold

    main:
    ch_multiqc_files  = Channel.empty()
    ch_top_ranked_pdb = Channel.empty()
    ch_pdb            = Channel.empty()
    ch_plddt          = Channel.empty()
    ch_msa            = Channel.empty()
    ch_paes           = Channel.empty()
    ch_multiqc_report = Channel.empty()

    //
    // MODULE: Run esmfold
    //
    if (params.esmfold_model_preset != 'monomer') {
        MULTIFASTA_TO_SINGLEFASTA(
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(MULTIFASTA_TO_SINGLEFASTA.out.versions)
        RUN_ESMFOLD(
            MULTIFASTA_TO_SINGLEFASTA.out.input_fasta,
            ch_esmfold_params,
            ch_num_recycles
        )

        ch_top_ranked_pdb = ch_top_ranked_pdb.mix(RUN_ESMFOLD.out.top_ranked_pdb)
        ch_pdb            = ch_pdb.mix(RUN_ESMFOLD.out.pdb)
        ch_plddt          = ch_plddt.mix(RUN_ESMFOLD.out.plddt)
        ch_msa            = ch_msa.mix(RUN_ESMFOLD.out.msa)
        ch_paes           = ch_paes.mix(RUN_ESMFOLD.out.paes)
        ch_versions       = ch_versions.mix(RUN_ESMFOLD.out.versions)

    } else {
        RUN_ESMFOLD(
            ch_samplesheet,
            ch_esmfold_params,
            ch_num_recycles
        )

        ch_top_ranked_pdb = ch_top_ranked_pdb.mix(RUN_ESMFOLD.out.top_ranked_pdb)
        ch_pdb            = ch_pdb.mix(RUN_ESMFOLD.out.pdb)
        ch_plddt          = ch_plddt.mix(RUN_ESMFOLD.out.plddt)
        ch_msa            = ch_msa.mix(RUN_ESMFOLD.out.msa)
        ch_paes           = ch_paes.mix(RUN_ESMFOLD.out.paes)
        ch_versions       = ch_versions.mix(RUN_ESMFOLD.out.versions)

    }

    emit:
    top_ranked_pdb = ch_top_ranked_pdb // channel: [ id, /path/to/*.pdb ]
    pdb            = ch_pdb
    plddt          = ch_plddt
    msa            = ch_msa
    paes           = ch_paes
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions         // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
