#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteinfold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteinfold
    Website: https://nf-co.re/proteinfold
    Slack  : https://nfcore.slack.com/channels/proteinfold
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_ALPHAFOLD2_DBS           } from './subworkflows/local/prepare_alphafold2_dbs'
include { PREPARE_ALPHAFOLD3_DBS           } from './subworkflows/local/prepare_alphafold3_dbs'

include { PREPARE_COLABFOLD2_DBS           } from './subworkflows/local/prepare_colabfold2_dbs'

include { PROTEINFOLD                      } from './workflows/proteinfold'

include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { modelcifInput                    } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { POST_PROCESSING                  } from './subworkflows/local/post_processing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//

workflow NFCORE_PROTEINFOLD {

    take:
    samplesheet  // channel: samplesheet read in from --input

    main:
    ch_samplesheet       = samplesheet
    ch_multiqc           = channel.empty()
    ch_report_input      = channel.empty()
    ch_top_ranked_model  = channel.empty()
    modelcif_software_details = file("$projectDir/assets/DUMMY_SOFTWARE_DETAILS.yml", checkIfExists: true)
    def ch_software_versions = channel.topic('versions')
        .unique()
        .map { process_name, tool_name, version ->
            "\"${process_name}:${tool_name}\": ${version}"
        }
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_proteinfold_software_mqc_versions.yml',
            newLine: true,
            sort: true
        )
    requested_modes      = params.mode.toLowerCase().split(",")

    colabfold2_prefixed_models = [ 'alphafold3', 'boltz2' ]
    colabfold2_bare_models     = [
        'openfold3', 'openbind0', 'protenix2', 'chai1', 'intellifold2', 'opendde',
        'rosettafold3', 'esmfold2',
        'alphafold2'
    ]
    colabfold2_mode_tokens = ( colabfold2_prefixed_models.collect { "colabfold2-${it}".toString() } + colabfold2_bare_models )
    colabfold2_submodes = requested_modes
        .collect { mode -> mode.trim() }
        .findAll { mode -> colabfold2_mode_tokens.contains(mode) }
    if (requested_modes.collect { it.trim() }.contains('colabfold2')) {
        error("The bare `colabfold2` mode selects no backend. Name one, e.g. `--mode colabfold2-alphafold3`, or one of: ${colabfold2_bare_models.join(', ')}.")
    }
    colabfold2_lookalikes = requested_modes
        .collect { mode -> mode.trim() }
        .findAll { mode -> mode.startsWith('colabfold2') || mode.startsWith('colabfold-') }
    colabfold2_unknown_submodes = requested_modes
        .collect { mode -> mode.trim() }
        .findAll { mode ->
            mode != 'colabfold2' &&
            !colabfold2_submodes.contains(mode) &&
            ( colabfold2_lookalikes.contains(mode) || ( colabfold2_bare_models + colabfold2_prefixed_models ).any { backend -> mode.startsWith(backend) && mode != backend } )
        }
    if (colabfold2_unknown_submodes) {
        error("Unknown ColabFold2 mode(s): ${colabfold2_unknown_submodes.join(', ')}. Supported modes: ${colabfold2_mode_tokens.join(', ')}.")
    }
    if (colabfold2_submodes.size() > 1) {
        error("Running multiple ColabFold2 backends in one invocation is not supported yet; run them separately.")
    }
    colabfold2_public_mode = colabfold2_submodes ? colabfold2_submodes[0] : null
    colabfold2_esmfold2_models = [
        'msa'    : [ model_type: 'esmfold2',        use_esm: false ],
        'lm300m' : [ model_type: 'esmfold2_lm300m', use_esm: true  ],
        'lm600m' : [ model_type: 'esmfold2_lm600m', use_esm: true  ]
    ]
    colabfold2_esmfold2_selection = params.esmfold2_model ?: 'msa'
    if (colabfold2_public_mode == 'colabfold2-alphafold3' && params.colabfold2_weights_precision == 'int8') {
        error("ColabFold2 official alphafold3 supports fp32 weights only; --colabfold2_weights_precision int8 cannot be used with --mode colabfold2-alphafold3.")
    }
    if (colabfold2_public_mode && colabfold2_public_mode.startsWith('alphafold2') && params.colabfold2_weights_precision == 'int8') {
        error("The AlphaFold2 backends do not provide int8 weights; --colabfold2_weights_precision int8 cannot be used with --mode ${colabfold2_public_mode}.")
    }
    if (colabfold2_public_mode == 'esmfold2') {
        def esmfold2_model = colabfold2_esmfold2_models[colabfold2_esmfold2_selection]
        if (esmfold2_model == null) {
            error("Unsupported --esmfold2_model '${colabfold2_esmfold2_selection}'. Choose one of: ${colabfold2_esmfold2_models.keySet().join(', ')}.")
        }
        colabfold2_model_for_submode = esmfold2_model.model_type
        colabfold2_use_esm = esmfold2_model.use_esm
    } else {
        colabfold2_model_for_submode = colabfold2_public_mode ?
            colabfold2_public_mode.replaceFirst(/^colabfold2-/, '') : null
        colabfold2_use_esm = false
        // The plain alphafold2 token delegates monomer/multimer selection to
        // colabfold_batch ('auto': alphafold2_ptm / alphafold2_multimer_v3).
        if (colabfold2_model_for_submode == 'alphafold2') {
            colabfold2_model_for_submode = 'auto'
        }
    }
    if (colabfold2_submodes && !(colabfold2_model_for_submode in ( colabfold2_prefixed_models + colabfold2_bare_models + [ 'esmfold2_lm300m', 'esmfold2_lm600m', 'auto' ] ))) {
        error("ColabFold2 mode '${colabfold2_submodes[0]}' resolved to backend '${colabfold2_model_for_submode}', which is not a supported AlphaFold3-family backend. Supported modes: ${colabfold2_mode_tokens.join(', ')}.")
    }

    //
    // Resolve the MSA generation module
    //
    msa_tool = params.use_msa_server ? 'msa_server' : params.msa_tool.toLowerCase()
    msa_tools = ['colabfold_search', 'alphafold2_msa', 'alphafold3_msa']
    if (!params.use_msa_server && !msa_tools.contains(msa_tool)) {
        error("Unknown --msa_tool '${msa_tool}'. Choose one of: ${msa_tools.join(', ')} or set --use_msa_server.")
    }

    ch_af2_params       = channel.empty()
    ch_af2_bfd          = channel.empty()
    ch_af2_small_bfd    = channel.empty()
    ch_af2_mgnify       = channel.empty()
    ch_af2_pdb70        = channel.empty()
    ch_af2_pdb_mmcif    = channel.empty()
    ch_af2_pdb_obsolete = channel.empty()
    ch_af2_uniref30     = channel.empty()
    ch_af2_uniref90     = channel.empty()
    ch_af2_pdb_seqres   = channel.empty()
    ch_af2_uniprot      = channel.empty()

    ch_af3_small_bfd    = channel.empty()
    ch_af3_mgnify       = channel.empty()
    ch_af3_mmcif_files  = channel.empty()
    ch_af3_uniref90     = channel.empty()
    ch_af3_pdb_seqres   = channel.empty()
    ch_af3_uniprot      = channel.empty()
    ch_af3_nt_rna       = channel.empty()
    ch_af3_rfam         = channel.empty()
    ch_af3_rnacentral   = channel.empty()

    if (msa_tool == 'alphafold2_msa') {

        //
        // SUBWORKFLOW: Prepare Alphafold2 DBs
        //
        PREPARE_ALPHAFOLD2_DBS (
            params.alphafold2_db,
            params.alphafold2_full_dbs,
            params.alphafold2_bfd_path,
            params.alphafold2_small_bfd_path,
            params.alphafold2_params_path,
            params.alphafold2_mgnify_path,
            params.alphafold2_pdb70_path,
            params.alphafold2_pdb_mmcif_path,
            params.alphafold2_pdb_obsolete_path,
            params.alphafold2_uniref30_path,
            params.alphafold2_uniref90_path,
            params.alphafold2_pdb_seqres_path,
            params.alphafold2_uniprot_path,
            params.alphafold2_bfd_link,
            params.alphafold2_small_bfd_link,
            params.alphafold2_params_link,
            params.alphafold2_mgnify_link,
            params.alphafold2_pdb70_link,
            params.alphafold2_pdb_mmcif_link,
            params.alphafold2_pdb_obsolete_link,
            params.alphafold2_uniref30_link,
            params.alphafold2_uniref90_link,
            params.alphafold2_pdb_seqres_link,
            params.alphafold2_uniprot_sprot_link,
            params.alphafold2_uniprot_trembl_link
        )

        ch_af2_params       = PREPARE_ALPHAFOLD2_DBS.out.params
        ch_af2_bfd          = PREPARE_ALPHAFOLD2_DBS.out.bfd
        ch_af2_small_bfd    = PREPARE_ALPHAFOLD2_DBS.out.small_bfd
        ch_af2_mgnify       = PREPARE_ALPHAFOLD2_DBS.out.mgnify
        ch_af2_pdb70        = PREPARE_ALPHAFOLD2_DBS.out.pdb70
        ch_af2_pdb_mmcif    = PREPARE_ALPHAFOLD2_DBS.out.pdb_mmcif
        ch_af2_pdb_obsolete = PREPARE_ALPHAFOLD2_DBS.out.pdb_obsolete
        ch_af2_uniref30     = PREPARE_ALPHAFOLD2_DBS.out.uniref30
        ch_af2_uniref90     = PREPARE_ALPHAFOLD2_DBS.out.uniref90
        ch_af2_pdb_seqres   = PREPARE_ALPHAFOLD2_DBS.out.pdb_seqres
        ch_af2_uniprot      = PREPARE_ALPHAFOLD2_DBS.out.uniprot
    }

    if (msa_tool == 'alphafold3_msa') {

        //
        // SUBWORKFLOW: Prepare Alphafold3 DBs
        //
        PREPARE_ALPHAFOLD3_DBS (
            params.alphafold3_db,
            params.alphafold3_params_path,
            params.alphafold3_small_bfd_path,
            params.alphafold3_mgnify_path,
            params.alphafold3_pdb_mmcif_path,
            params.alphafold3_uniref90_path,
            params.alphafold3_pdb_seqres_path,
            params.alphafold3_uniprot_path,
            params.alphafold3_rnacentral_path,
            params.alphafold3_nt_rna_path,
            params.alphafold3_rfam_path,
            params.alphafold3_small_bfd_link,
            params.alphafold3_mgnify_link,
            params.alphafold3_pdb_mmcif_link,
            params.alphafold3_uniref90_link,
            params.alphafold3_pdb_seqres_link,
            params.alphafold3_uniprot_link,
            params.alphafold3_rnacentral_link,
            params.alphafold3_nt_rna_link,
            params.alphafold3_rfam_link
        )

        ch_af3_small_bfd    = PREPARE_ALPHAFOLD3_DBS.out.small_bfd
        ch_af3_mgnify       = PREPARE_ALPHAFOLD3_DBS.out.mgnify
        ch_af3_mmcif_files  = PREPARE_ALPHAFOLD3_DBS.out.pdb_mmcif
        ch_af3_uniref90     = PREPARE_ALPHAFOLD3_DBS.out.uniref90
        ch_af3_pdb_seqres   = PREPARE_ALPHAFOLD3_DBS.out.pdb_seqres
        ch_af3_uniprot      = PREPARE_ALPHAFOLD3_DBS.out.uniprot
        ch_af3_nt_rna       = PREPARE_ALPHAFOLD3_DBS.out.nt_rna
        ch_af3_rfam         = PREPARE_ALPHAFOLD3_DBS.out.rfam
        ch_af3_rnacentral   = PREPARE_ALPHAFOLD3_DBS.out.rnacentral
    }

    if (!params.colabfold2_params_path && !params.colabfold_db && !params.use_msa_server) {
        error("ColabFold2 needs weights: provide --colabfold2_params_path, --db/--colabfold_db, or use --use_msa_server. Weights are downloaded separately when no local weights root is provided.")
    }

    //
    // SUBWORKFLOW: Prepare ColabFold2 weights/DBs
    //
    PREPARE_COLABFOLD2_DBS (
        params.colabfold2_params_path,
        params.colabfold_db,
        colabfold2_model_for_submode,
        params.colabfold2_weights_precision,
        params.use_msa_server,
        params.colabfold_envdb_path,
        params.colabfold_uniref30_path,
        params.colabfold_db_link,
        params.colabfold_uniref30_link,
        params.colabfold_create_index
    )
    ch_colabfold2_params_path = PREPARE_COLABFOLD2_DBS.out.params_path
    ch_colabfold2_db          = PREPARE_COLABFOLD2_DBS.out.colabfold_db
    ch_colabfold2_uniref30    = PREPARE_COLABFOLD2_DBS.out.uniref30

    //
    // WORKFLOW: Run the unified proteinfold workflow
    //
    PROTEINFOLD (
        ch_samplesheet,
        ch_colabfold2_params_path,
        colabfold2_model_for_submode,
        colabfold2_use_esm,
        colabfold2_public_mode,
        params.colabfold_num_recycles,
        msa_tool,
        params.msa_batch_size,
        ch_colabfold2_db,
        ch_colabfold2_uniref30,
        params.alphafold2_full_dbs != null ? params.alphafold2_full_dbs : false,
        params.uniref30_prefix ?: 'UniRef30_2023_02',
        ch_af2_params,
        ch_af2_bfd,
        ch_af2_small_bfd,
        ch_af2_mgnify,
        ch_af2_pdb70,
        ch_af2_pdb_mmcif,
        ch_af2_pdb_obsolete,
        ch_af2_uniref30,
        ch_af2_uniref90,
        ch_af2_pdb_seqres,
        ch_af2_uniprot,
        ch_af3_small_bfd,
        ch_af3_mgnify,
        ch_af3_mmcif_files,
        ch_af3_uniref90,
        ch_af3_pdb_seqres,
        ch_af3_uniprot,
        ch_af3_nt_rna,
        ch_af3_rfam,
        ch_af3_rnacentral
    )

    ch_multiqc      = ch_multiqc.mix(PROTEINFOLD.out.multiqc_metrics)
    ch_report_input = ch_report_input.mix(
        PROTEINFOLD.out.pdb
            .join(PROTEINFOLD.out.msa)
            .join(PROTEINFOLD.out.pae)
            .join(PROTEINFOLD.out.iptm)
            .join(PROTEINFOLD.out.ipsae)
            .join(PROTEINFOLD.out.chainwise_iptm)
            .join(PROTEINFOLD.out.chainwise_ipsae)
    )
    ch_top_ranked_model = ch_top_ranked_model.mix(PROTEINFOLD.out.top_ranked_pdb)

    //
    // POST PROCESSING: generate visualisation reports
    //
    ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ).first()  : channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo ).first()    : channel.empty()
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_report_template     = channel.value(file("$projectDir/assets/report_template.html", checkIfExists: true))

    // Inject msa_tool into meta based on selected MSA module for provenance
    // (report and ModelCIF assembly).
    def msaToolMap = [
        colabfold_search: 'mmseqs2',
        alphafold2_msa:   'jackhmmer',
        alphafold3_msa:   'jackhmmer',
        msa_server:       'mmseqs2-server',
    ]
    def selected_msa_tool = msaToolMap.get(msa_tool, 'None')
    ch_report_input = ch_report_input.map { tupleData ->
        def meta = tupleData[0]
        def m = meta.clone()
        m.msa_tool = selected_msa_tool
        [m] + tupleData.drop(1)
    }

    ch_modelcif = modelcifInput(
        PROTEINFOLD.out.pdb,
        PROTEINFOLD.out.msa,
        PROTEINFOLD.out.plddt,
        PROTEINFOLD.out.pae,
        PROTEINFOLD.out.ptm,
        PROTEINFOLD.out.iptm,
        ch_software_versions,
        modelcif_software_details,
        selected_msa_tool
    )

    POST_PROCESSING(
        params.skip_visualisation,
        ch_report_input,
        ch_report_template,
        params.skip_foldseek,
        params.foldseek_db,
        params.foldseek_db_path,
        params.skip_multiqc,
        params.outdir,
        ch_multiqc,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        params.multiqc_logo,
        ch_multiqc_methods_description,
        ch_software_versions,
        ch_top_ranked_model,
        ch_modelcif,
        params.write_modelcif,
        params.modelcif_binary,
        params.modelcif_pae_embed,
        params.modelcif_software_details
    )

    emit:
    multiqc_report = POST_PROCESSING.out.multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PROTEINFOLD (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_PROTEINFOLD.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
