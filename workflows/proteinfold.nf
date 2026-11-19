/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified proteinfold workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Single workflow for the pipeline: users choose one MSA generation module
    (--msa_tool: colabfold_search | alphafold2_msa | alphafold3_msa, or
    --use_msa_server for the ColabFold webserver), and every structure
    prediction step runs through the colabfold2 backends.
*/

include { COLABFOLD2_INFERENCE as ALPHAFOLD2      } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as ALPHAFOLD3      } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as BOLTZ2          } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as OPENFOLD3       } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as OPENBIND0       } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as PROTEINX2       } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as CHAI1           } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as INTELLIFOLD2    } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as OPENDDE         } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as ROSETTAFOLD3    } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE as ESMFOLD2        } from '../modules/local/colabfold2_inference'

include { ALPHAFOLD3_JSON_TO_COLABFOLD_FASTA } from '../modules/local/alphafold3_json_to_colabfold_fasta'
include { ALPHAFOLD3_JSON_TO_COLABFOLD_A3M   } from '../modules/local/af3_msa_json_to_a3m'
include { MSA_COLABFOLD_SEARCH } from '../subworkflows/local/msa_colabfold_search'
include { MSA_ALPHAFOLD2       } from '../subworkflows/local/msa_alphafold2'
include { MSA_ALPHAFOLD3       } from '../subworkflows/local/msa_alphafold3'
include { collectMultiqcMetrics; modeChannel } from '../subworkflows/local/utils_nfcore_proteinfold_pipeline'

workflow PROTEINFOLD {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    ch_params_path  // channel: path(colabfold2_params_path, the weights/data dir root)
    model_type      // string: explicit upstream AF3-family model name
    use_esm         // boolean: select the ESM-C language-model path for ESMFold2
    public_mode     // string: public mode token used for output publication
    num_recycles    // int: Number of recycles for the selected backend
    msa_tool        // string: colabfold_search | alphafold2_msa | alphafold3_msa | msa_server
    msa_batch_size  // int: number of samples processed together during MSA search (-1 = one batch)
    ch_colabfold_db // channel: path(colabfold_db)
    ch_uniref30     // channel: path(uniref30)
    // alphafold2_msa tool inputs
    alphafold2_full_dbs // boolean: full AF2 databases (otherwise reduced)
    alphafold2_max_template_date // string: latest template release date to keep
    uniref30_prefix     // string: prefix for uniref30 database files
    ch_af2_params       // channel: path(alphafold2 params)
    ch_af2_bfd          // channel: path(bfd)
    ch_af2_small_bfd    // channel: path(small_bfd)
    ch_af2_mgnify       // channel: path(mgnify)
    ch_af2_pdb70        // channel: path(pdb70)
    ch_af2_pdb_mmcif    // channel: path(pdb_mmcif)
    ch_af2_pdb_obsolete // channel: path(pdb_obsolete)
    ch_af2_uniref30     // channel: path(uniref30)
    ch_af2_uniref90     // channel: path(uniref90)
    ch_af2_pdb_seqres   // channel: path(pdb_seqres)
    ch_af2_uniprot      // channel: path(uniprot)
    // alphafold3_msa tool inputs
    ch_af3_small_bfd    // channel: path(small_bfd)
    ch_af3_mgnify       // channel: path(mgnify)
    ch_af3_mmcif_files  // channel: path(mmcif_files)
    ch_af3_uniref90     // channel: path(uniref90)
    ch_af3_pdb_seqres   // channel: path(pdb_seqres)
    ch_af3_uniprot      // channel: path(uniprot)
    ch_af3_nt_rna       // channel: path(nt_rna)
    ch_af3_rfam         // channel: path(rfam)
    ch_af3_rnacentral   // channel: path(rnacentral)

    main:
    def COLABFOLD2_BACKENDS = [
        'auto'                   : 'ALPHAFOLD2',   // the plain alphafold2 token
        'alphafold2'             : 'ALPHAFOLD2',   // defensive alias for the same token
        'alphafold3'      : 'ALPHAFOLD3',
        'boltz2'          : 'BOLTZ2',
        'openfold3'       : 'OPENFOLD3',
        'openbind0'       : 'OPENBIND0',
        'protenix2'       : 'PROTEINX2',
        'chai1'           : 'CHAI1',
        'intellifold2'    : 'INTELLIFOLD2',
        'opendde'         : 'OPENDDE',
        'rosettafold3'    : 'ROSETTAFOLD3',
        'esmfold2'        : 'ESMFOLD2',
        'esmfold2_lm300m' : 'ESMFOLD2',
        'esmfold2_lm600m' : 'ESMFOLD2'
    ]

    def backend_process = COLABFOLD2_BACKENDS[model_type]
    if (backend_process == null) {
        error("ColabFold2 backend '${model_type}' is not wired to an inference process. Wired backends: ${COLABFOLD2_BACKENDS.keySet().join(', ')}.")
    }

    ch_byo_msa = ch_samplesheet
        .filter { meta, input -> meta.msa }
        .map { meta, input ->
            def clean_meta = meta.clone()
            def msa = clean_meta.remove('msa')
            [clean_meta, msa]
        }

    ch_samplesheet
        .filter { meta, input -> !meta.msa }
        .map { meta, input ->
            def extension = input.extension.toLowerCase()
            if (!(extension in ['fa', 'faa', 'fasta', 'json'])) {
                error("Proteinfold input for sample '${meta.id}' must be FASTA or AlphaFold3 JSON; received '${input.name}'.")
            }
            [meta, input]
        }
        .branch { meta, input ->
            json:  input.extension.toLowerCase() == 'json'
            fasta: true
        }
        .set { ch_input_by_ext }

    //
    // MSA generation: one of the three MSA modules (or the ColabFold server)
    //

    if (msa_tool == 'msa_server') {
        // colabfold_batch performs the MSA search against the webserver itself.
        ch_msa_input = ch_input_by_ext.fasta.mix(ch_input_by_ext.json)
    } else if (msa_tool == 'colabfold_search') {
        MSA_COLABFOLD_SEARCH(
            ch_input_by_ext.fasta,
            ch_input_by_ext.json,
            msa_batch_size,
            ch_colabfold_db,
            ch_uniref30
        )
        ch_msa_input = MSA_COLABFOLD_SEARCH.out.msa_json
    } else if (msa_tool == 'alphafold2_msa') {
        MSA_ALPHAFOLD2(
            ch_input_by_ext.fasta,
            ch_input_by_ext.json,
            alphafold2_full_dbs,
            alphafold2_max_template_date,
            uniref30_prefix,
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
            ch_af2_uniprot
        )
        ch_msa_input = MSA_ALPHAFOLD2.out.msa_json
    } else if (msa_tool == 'alphafold3_msa') {
        MSA_ALPHAFOLD3(
            ch_input_by_ext.fasta,
            ch_input_by_ext.json,
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
        ch_msa_input = MSA_ALPHAFOLD3.out.msa_json
    } else {
        error("Unknown --msa_tool '${msa_tool}'. Choose one of: colabfold_search, alphafold2_msa, alphafold3_msa, or set --use_msa_server.")
    }

    //
    // Structure prediction: colabfold2 backends
    //

    // The alphafold2 backends predate the AlphaFold3 input JSON: they consume
    // ColabFold a3m queries (or a plain FASTA when the MSA server does the
    // search). Convert whichever MSA JSON the selected tool produced.
    def af2_backend = model_type == 'auto' || model_type.startsWith('alphafold2') || model_type.startsWith('deepfold')
    if (af2_backend && msa_tool != 'msa_server') {
        ALPHAFOLD3_JSON_TO_COLABFOLD_A3M(ch_msa_input)
        ch_inference_input = ALPHAFOLD3_JSON_TO_COLABFOLD_A3M.out.a3m.mix(ch_byo_msa)
    } else if (af2_backend) {
        // Server mode: let colabfold_batch run the server search itself; JSON
        // inputs are reduced to their protein chains first.
        ALPHAFOLD3_JSON_TO_COLABFOLD_FASTA(ch_input_by_ext.json)
        ch_inference_input = ch_input_by_ext.fasta.mix(ALPHAFOLD3_JSON_TO_COLABFOLD_FASTA.out.query_fasta).mix(ch_byo_msa)
    } else {
        ch_inference_input = ch_msa_input.mix(ch_byo_msa)
    }

    if (backend_process == 'ALPHAFOLD2') {
        ALPHAFOLD2(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = ALPHAFOLD2.out
        ch_batch_structures = ALPHAFOLD2.out.structures_pdb
    } else if (backend_process == 'ALPHAFOLD3') {
        ALPHAFOLD3(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = ALPHAFOLD3.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'BOLTZ2') {
        BOLTZ2(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = BOLTZ2.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'OPENFOLD3') {
        OPENFOLD3(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = OPENFOLD3.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'OPENBIND0') {
        OPENBIND0(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = OPENBIND0.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'PROTEINX2') {
        PROTEINX2(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = PROTEINX2.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'CHAI1') {
        CHAI1(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = CHAI1.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'INTELLIFOLD2') {
        INTELLIFOLD2(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = INTELLIFOLD2.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'OPENDDE') {
        OPENDDE(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = OPENDDE.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'ROSETTAFOLD3') {
        ROSETTAFOLD3(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = ROSETTAFOLD3.out
        ch_batch_structures = ch_batch.structures
    } else if (backend_process == 'ESMFOLD2') {
        ESMFOLD2(ch_inference_input, model_type, use_esm, public_mode, num_recycles, ch_params_path)
        ch_batch = ESMFOLD2.out
        ch_batch_structures = ch_batch.structures
    } else {
        error("ColabFold2 backend '${model_type}' resolved to unknown process '${backend_process}'.")
    }

    ch_batch.top_ranked_structure
        .map { meta, structure ->
            def report_meta = meta.clone()
            report_meta.model = 'colabfold2'
            [report_meta, structure]
        }
        .set { ch_top_ranked_structure }

    ch_batch_structures
        .map { meta, structures ->
            def report_meta = meta.clone()
            report_meta.model = 'colabfold2'
            def ranked = (structures instanceof List ? structures : [structures]).sort { structure ->
                def matcher = structure.name =~ /_rank_(\d+)_/
                matcher.find() ? matcher.group(1).toInteger() : Integer.MAX_VALUE
            }
            [report_meta, ranked]
        }
        .set { ch_structures }

    modeChannel(ch_batch.msa, 'colabfold2').set { ch_msa }
    modeChannel(ch_batch.pae, 'colabfold2').set { ch_pae }
    modeChannel(ch_batch.ptms, 'colabfold2').set { ch_ptm }
    modeChannel(ch_batch.iptms, 'colabfold2').set { ch_iptm }
    modeChannel(ch_batch.ipsaes, 'colabfold2').set { ch_ipsae }
    modeChannel(ch_batch.chainwise_iptms, 'colabfold2').set { ch_chainwise_iptm }
    modeChannel(ch_batch.chainwise_ipsaes, 'colabfold2').set { ch_chainwise_ipsae }
    modeChannel(ch_batch.confidence, 'colabfold2').set { ch_confidence }
    modeChannel(ch_batch.plddt, 'colabfold2').set { ch_plddt }

    ch_multiqc_metrics = collectMultiqcMetrics('colabfold2', [
        [ 'plddt', ch_batch.plddt ],
        [ 'msa',   ch_batch.msa ],
        [ 'ptm',   ch_batch.ptms ],
        [ 'iptm',  ch_batch.iptms ]
    ])

    emit:
    top_ranked_pdb  = ch_top_ranked_structure
    pdb             = ch_structures
    msa             = ch_msa
    pae             = ch_pae
    ptm             = ch_ptm
    iptm            = ch_iptm
    ipsae           = ch_ipsae
    chainwise_iptm  = ch_chainwise_iptm
    chainwise_ipsae = ch_chainwise_ipsae
    confidence      = ch_confidence
    plddt           = ch_plddt
    multiqc_metrics = ch_multiqc_metrics
}
