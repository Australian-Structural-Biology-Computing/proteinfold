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
include { COLABFOLD2_INFERENCE_BATCH as ALPHAFOLD2_BATCH } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as ALPHAFOLD3_BATCH } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as BOLTZ2_BATCH     } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as OPENFOLD3_BATCH  } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as OPENBIND0_BATCH  } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as PROTEINX2_BATCH  } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as CHAI1_BATCH      } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as INTELLIFOLD2_BATCH } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as OPENDDE_BATCH    } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as ROSETTAFOLD3_BATCH } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as ESMFOLD2_BATCH   } from '../modules/local/colabfold2_inference'

include { TO_COLABFOLD_FASTA } from '../modules/local/to_colabfold_fasta'
include { ALPHAFOLD3_JSON_TO_COLABFOLD_A3M   } from '../modules/local/af3_msa_json_to_a3m'
include { MSA_COLABFOLD_SEARCH } from '../subworkflows/local/msa_colabfold_search'
include { MSA_ALPHAFOLD2       } from '../subworkflows/local/msa_alphafold2'
include { MSA_ALPHAFOLD3       } from '../subworkflows/local/msa_alphafold3'
include { collectMultiqcMetrics; modeChannel } from '../subworkflows/local/utils_nfcore_proteinfold_pipeline'


// Map a file produced by COLABFOLD2_INFERENCE_BATCH to its sample by matching
// the longest sample-id prefix (sample ids may be prefixes of one another).
def matchBatchSample(batch_metas, filename) {
    def ids = batch_metas.collect { it.id as String }.sort { -it.length() }
    def id = ids.find { filename.startsWith(it + '_') || filename.startsWith(it + '.') }
    if (id == null) {
        error("COLABFOLD2_INFERENCE_BATCH produced file '${filename}' that does not match any sample in the batch (${ids.join(', ')})")
    }
    return batch_metas.find { (it.id as String) == id }
}

// One file per sample: emit [meta, file] tuples.
def splitBatchSingleFiles(out_channel) {
    out_channel.flatMap { batch_meta, files ->
        def list = files instanceof List ? files : [files]
        list.collect { file -> [matchBatchSample(batch_meta, file.name), file] }
    }
}

// Several files per sample (ranked structures, score JSONs): emit [meta, [files]] sorted by rank.
def splitBatchGroupedFiles(out_channel) {
    out_channel.flatMap { batch_meta, files ->
        def list = files instanceof List ? files : [files]
        def grouped = list
            .collect { file -> [matchBatchSample(batch_meta, file.name), file] }
            .groupBy { it[0] }
        batch_meta.collect { meta ->
            def sample_files = (grouped[meta] ?: []).collect { it[1] }.sort { file ->
                def matcher = file.name =~ /_rank_(\d+)_/
                matcher.find() ? matcher.group(1).toInteger() : Integer.MAX_VALUE
            }
            [meta, sample_files]
        }
    }
}

// Normalise the batch process outputs into the per-sample channel map the
// single-sample dispatch produces.
def splitBatchOutputs(out) {
    [
        top_ranked_structure: splitBatchSingleFiles(out.top_ranked_structures),
        structures          : splitBatchGroupedFiles(out.structures_cif),
        structures_pdb      : splitBatchGroupedFiles(out.structures_pdb),
        msa                 : splitBatchSingleFiles(out.msa),
        pae                 : splitBatchSingleFiles(out.pae),
        ptms                : splitBatchSingleFiles(out.ptms),
        iptms               : splitBatchSingleFiles(out.iptms),
        ipsaes              : splitBatchSingleFiles(out.ipsaes),
        chainwise_iptms     : splitBatchSingleFiles(out.chainwise_iptms),
        chainwise_ipsaes    : splitBatchSingleFiles(out.chainwise_ipsaes),
        confidence          : splitBatchSingleFiles(out.confidence),
        plddt               : splitBatchSingleFiles(out.plddt),
        scores              : splitBatchGroupedFiles(out.scores),
        raw                 : out.raw
    ]
}

// Group one backend's per-sample stream into batch tasks. Samples are grouped
// by model type, ESM flag and weights directory (uniform per colabfold_batch
// invocation), sorted by id inside each group, and collated into batches of
// --inference_batch_size (-1 = one batch of everything).
def batchInferenceInputs(ch, inference_batch) {
    ch
        .map { meta, input, model_type, use_esm, dir ->
            tuple("${model_type}|${use_esm}|${dir}", meta, input, model_type, use_esm, dir)
        }
        .groupTuple()
        .flatMap { key, metas, inputs, model_types, use_esms, dirs ->
            def sorted = [metas, inputs].transpose().sort { a, b -> (a[0].id as String) <=> (b[0].id as String) }
            def effective = inference_batch == -1 ? Math.max(1, sorted.size()) : inference_batch
            sorted.collate(effective).collect { batch ->
                tuple(batch*.getAt(0) as List, batch*.getAt(1) as List, model_types[0], use_esms[0], dirs[0])
            }
        }
}

workflow PROTEINFOLD {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    ch_params_path  // channel: tuple(model_type, params dir) from PREPARE_COLABFOLD2_DBS
    mode_specs      // list of maps: [public_mode, model_type, use_esm] per requested --mode
    num_recycles    // int: Number of recycles for the selected backend
    msa_tool        // string: colabfold_search | alphafold2_msa | alphafold3_msa | msa_server
    msa_batch_size  // int: number of samples processed together during MSA search (-1 = one batch)
    inference_batch_size // int: number of samples processed together during structure prediction (-1 = one batch)
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
    // The compiled-kernel cache is bound into the inference containers, and
    // a bind source must exist before a container starts -- create it on the
    // host, owned by the launching user, ahead of any task.
    if (params.colabfold2_compile_cache) {
        // Resolved against the launch directory, matching the module's bind.
        new File(file(params.colabfold2_compile_cache).toAbsolutePath().toString()).mkdirs()
    }
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

    ch_byo_msa = ch_samplesheet
        .filter { meta, input -> meta.msa }
        .map { meta, input ->
            def clean_meta = meta.clone()
            def msa = clean_meta.remove('msa')
            clean_meta.msa_short = 'byo'
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
        ch_msa_input = ch_input_by_ext.fasta.mix(ch_input_by_ext.json).map { meta, input ->
            def tagged = meta.clone()
            tagged.msa_short = 'server'
            [tagged, input]
        }
    } else if (msa_tool == 'colabfold_search') {
        MSA_COLABFOLD_SEARCH(
            ch_input_by_ext.fasta,
            ch_input_by_ext.json,
            msa_batch_size,
            ch_colabfold_db,
            ch_uniref30
        )
        ch_msa_input = MSA_COLABFOLD_SEARCH.out.msa_json.map { meta, input ->
            def tagged = meta.clone()
            tagged.msa_short = 'mmseqs'
            [tagged, input]
        }
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
        ch_msa_input = MSA_ALPHAFOLD2.out.msa_json.map { meta, input ->
            def tagged = meta.clone()
            tagged.msa_short = 'af2msa'
            [tagged, input]
        }
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
        ch_msa_input = MSA_ALPHAFOLD3.out.msa_json.map { meta, input ->
            def tagged = meta.clone()
            tagged.msa_short = 'af3msa'
            [tagged, input]
        }
    } else {
        error("Unknown --msa_tool '${msa_tool}'. Choose one of: colabfold_search, alphafold2_msa, alphafold3_msa, or set --use_msa_server.")
    }

    //
    // Structure prediction: colabfold2 backends
    //

    // The alphafold2 backends predate the AlphaFold3 input JSON: they consume
    // ColabFold a3m queries (or a plain FASTA when the MSA server does the
    // search). Convert whichever MSA JSON the selected tool produced.

    // Split the requested modes into backend families. The AlphaFold2
    // backends run ColabFold's AF2 graph and consume ColabFold a3m queries:
    // they cannot fold non-protein entities. Every other mode runs the
    // AlphaFold3-open graph and accepts full AlphaFold3 JSON inputs,
    // including nucleic-acid and ligand entities (esmfold2 included).
    // Shorthand tokens baked into published file names, keeping the output
    // directory flat (no per-mode nesting).
    def MODE_SHORTS = [
        alphafold2: 'af2',
        alphafold3: 'af3',
        boltz2    : 'bz2',
        chai1     : 'chai1',
        esmfold2  : 'esm2',
        intellifold2: 'if2',
        opendde   : 'odde',
        openfold3 : 'of3',
        openbind0 : 'ob0',
        protenix2 : 'px2',
        rosettafold3: 'rf3',
    ]
    def af2_specs = mode_specs.findAll { it.model_type == 'auto' || it.model_type.startsWith('alphafold2') }
    def af3_specs = mode_specs.findAll { !(it.model_type == 'auto' || it.model_type.startsWith('alphafold2')) }

    def ch_af3_inference_input = null
    def ch_af2_inference_input = null

    if (af3_specs) {
        ch_af3_inference_input = ch_msa_input.mix(ch_byo_msa)
    }

    if (af2_specs) {
        if (msa_tool == 'msa_server') {
            // Server mode: colabfold_batch runs the webserver search itself;
            // JSON inputs are reduced to their protein chains first.
            TO_COLABFOLD_FASTA(ch_input_by_ext.json)
            ch_af2_inference_input = ch_input_by_ext.fasta
                .mix(TO_COLABFOLD_FASTA.out.query_fasta)
                .map { meta, input ->
                    def tagged = meta.clone()
                    tagged.msa_short = 'server'
                    [tagged, input]
                }
                .mix(ch_byo_msa)
        } else {
            // Convert the MSA JSON into ColabFold a3m queries. Samples whose
            // JSON cannot serve the AF2 graph (non-protein entities, no
            // protein entity, or no MSA to reuse) are flagged rather than
            // fatal so the remaining modes still fold them; the AF2 modes
            // skip those samples with a warning.
            ALPHAFOLD3_JSON_TO_COLABFOLD_A3M(ch_msa_input)
            def ch_convertible = ALPHAFOLD3_JSON_TO_COLABFOLD_A3M.out.a3m
                .map { meta, a3m, flag ->
                    def compatible = flag.name.startsWith('compatible')
                    if (!compatible) {
                        log.warn "Sample '${meta.id}' cannot be folded by the AlphaFold2 backends (${af2_specs*.public_mode.join(', ')}): its input carries non-protein entities or has no reusable protein MSA. It is still folded by: ${af3_specs ? af3_specs*.public_mode.join(', ') : 'no other requested mode'}."
                    }
                    [meta, a3m, compatible]
                }
            ch_af2_inference_input = ch_convertible
                .filter { meta, a3m, compatible -> compatible }
                .map { meta, a3m, compatible -> [meta, a3m] }
                .mix(ch_byo_msa)
        }
    }

    // Structure prediction: one invocation per backend process, at workflow
    // level (the Nextflow linter forbids process calls inside closures).
    // Each requested mode streams its samples into its backend's channel,
    // tagging the meta with the public mode and pairing every sample with
    // that mode's weights directory, so a single call covers every mode
    // routed to the backend.
    def ch_backends = [
        ALPHAFOLD2  : Channel.empty(),
        ALPHAFOLD3  : Channel.empty(),
        BOLTZ2      : Channel.empty(),
        OPENFOLD3   : Channel.empty(),
        OPENBIND0   : Channel.empty(),
        PROTEINX2   : Channel.empty(),
        CHAI1       : Channel.empty(),
        INTELLIFOLD2: Channel.empty(),
        OPENDDE     : Channel.empty(),
        ROSETTAFOLD3: Channel.empty(),
        ESMFOLD2    : Channel.empty()
    ]
    mode_specs.each { spec ->
        def backend = COLABFOLD2_BACKENDS[spec.model_type]
        if (backend == null) {
            error("ColabFold2 backend '${spec.model_type}' is not wired to an inference process. Wired backends: ${COLABFOLD2_BACKENDS.keySet().join(', ')}.")
        }
        def input_for_mode = (spec.model_type == 'auto' || spec.model_type.startsWith('alphafold2')) ? ch_af2_inference_input : ch_af3_inference_input
        // first() promotes the mode's weights directory to a value channel;
        // combine with a value channel pairs it with every sample (join
        // consumes right-side entries once, and combine against a queue is
        // a cartesian product).
        def params_for_mode = ch_params_path
            .filter { model, dir -> model == spec.model_type }
            .map { model, dir -> dir }
            .first()
        ch_backends[backend] = ch_backends[backend].mix(
            input_for_mode
                .combine(params_for_mode)
                .map { meta, input, dir ->
                    def tagged = meta.clone()
                    tagged.model = spec.public_mode
                    tagged.mode_short = MODE_SHORTS[spec.public_mode]
                    tuple(tagged, input, spec.model_type, spec.use_esm, dir)
                }
        )
    }

    // Inference batching mirrors --msa_batch_size: N samples share one
    // colabfold_batch task (one model/weight load for the whole batch).
    // -1 runs all samples of a backend's model in a single task; 1 (default)
    // runs per sample.
    def inference_batch = inference_batch_size instanceof Number ?
        inference_batch_size : inference_batch_size.toString().toInteger()
    if (inference_batch == 0 || inference_batch < -1) {
        error("Inference batch size must be -1 or a positive integer, got ${inference_batch}")
    }
    def batched_inference = inference_batch != 1
    if (batched_inference && msa_tool == 'msa_server') {
        log.warn "Inference batching is not supported with --use_msa_server; running inference per sample."
        batched_inference = false
    }

    def ch_backend_inputs = ch_backends.collectEntries { backend, ch ->
        [backend, batched_inference ? batchInferenceInputs(ch, inference_batch) : ch]
    }

    def backends_out = [:]
    if (af2_specs) {
        if (batched_inference) {
            ALPHAFOLD2_BATCH(ch_backend_inputs['ALPHAFOLD2'], num_recycles)
            backends_out['ALPHAFOLD2'] = splitBatchOutputs(ALPHAFOLD2_BATCH.out)
        } else {
            ALPHAFOLD2(ch_backend_inputs['ALPHAFOLD2'], num_recycles)
            backends_out['ALPHAFOLD2'] = ALPHAFOLD2.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'ALPHAFOLD3' }) {
        if (batched_inference) {
            ALPHAFOLD3_BATCH(ch_backend_inputs['ALPHAFOLD3'], num_recycles)
            backends_out['ALPHAFOLD3'] = splitBatchOutputs(ALPHAFOLD3_BATCH.out)
        } else {
            ALPHAFOLD3(ch_backend_inputs['ALPHAFOLD3'], num_recycles)
            backends_out['ALPHAFOLD3'] = ALPHAFOLD3.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'BOLTZ2' }) {
        if (batched_inference) {
            BOLTZ2_BATCH(ch_backend_inputs['BOLTZ2'], num_recycles)
            backends_out['BOLTZ2'] = splitBatchOutputs(BOLTZ2_BATCH.out)
        } else {
            BOLTZ2(ch_backend_inputs['BOLTZ2'], num_recycles)
            backends_out['BOLTZ2'] = BOLTZ2.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'OPENFOLD3' }) {
        if (batched_inference) {
            OPENFOLD3_BATCH(ch_backend_inputs['OPENFOLD3'], num_recycles)
            backends_out['OPENFOLD3'] = splitBatchOutputs(OPENFOLD3_BATCH.out)
        } else {
            OPENFOLD3(ch_backend_inputs['OPENFOLD3'], num_recycles)
            backends_out['OPENFOLD3'] = OPENFOLD3.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'OPENBIND0' }) {
        if (batched_inference) {
            OPENBIND0_BATCH(ch_backend_inputs['OPENBIND0'], num_recycles)
            backends_out['OPENBIND0'] = splitBatchOutputs(OPENBIND0_BATCH.out)
        } else {
            OPENBIND0(ch_backend_inputs['OPENBIND0'], num_recycles)
            backends_out['OPENBIND0'] = OPENBIND0.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'PROTEINX2' }) {
        if (batched_inference) {
            PROTEINX2_BATCH(ch_backend_inputs['PROTEINX2'], num_recycles)
            backends_out['PROTEINX2'] = splitBatchOutputs(PROTEINX2_BATCH.out)
        } else {
            PROTEINX2(ch_backend_inputs['PROTEINX2'], num_recycles)
            backends_out['PROTEINX2'] = PROTEINX2.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'CHAI1' }) {
        if (batched_inference) {
            CHAI1_BATCH(ch_backend_inputs['CHAI1'], num_recycles)
            backends_out['CHAI1'] = splitBatchOutputs(CHAI1_BATCH.out)
        } else {
            CHAI1(ch_backend_inputs['CHAI1'], num_recycles)
            backends_out['CHAI1'] = CHAI1.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'INTELLIFOLD2' }) {
        if (batched_inference) {
            INTELLIFOLD2_BATCH(ch_backend_inputs['INTELLIFOLD2'], num_recycles)
            backends_out['INTELLIFOLD2'] = splitBatchOutputs(INTELLIFOLD2_BATCH.out)
        } else {
            INTELLIFOLD2(ch_backend_inputs['INTELLIFOLD2'], num_recycles)
            backends_out['INTELLIFOLD2'] = INTELLIFOLD2.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'OPENDDE' }) {
        if (batched_inference) {
            OPENDDE_BATCH(ch_backend_inputs['OPENDDE'], num_recycles)
            backends_out['OPENDDE'] = splitBatchOutputs(OPENDDE_BATCH.out)
        } else {
            OPENDDE(ch_backend_inputs['OPENDDE'], num_recycles)
            backends_out['OPENDDE'] = OPENDDE.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'ROSETTAFOLD3' }) {
        if (batched_inference) {
            ROSETTAFOLD3_BATCH(ch_backend_inputs['ROSETTAFOLD3'], num_recycles)
            backends_out['ROSETTAFOLD3'] = splitBatchOutputs(ROSETTAFOLD3_BATCH.out)
        } else {
            ROSETTAFOLD3(ch_backend_inputs['ROSETTAFOLD3'], num_recycles)
            backends_out['ROSETTAFOLD3'] = ROSETTAFOLD3.out
        }
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'ESMFOLD2' }) {
        if (batched_inference) {
            ESMFOLD2_BATCH(ch_backend_inputs['ESMFOLD2'], num_recycles)
            backends_out['ESMFOLD2'] = splitBatchOutputs(ESMFOLD2_BATCH.out)
        } else {
            ESMFOLD2(ch_backend_inputs['ESMFOLD2'], num_recycles)
            backends_out['ESMFOLD2'] = ESMFOLD2.out
        }
    }

    // Combine every backend's emits. The metas already carry the public
    // mode (tagged at dispatch), so publication and MultiQC grouping stay
    // collision-free without further tagging.
    def combined = [top_ranked_structure: [], structures: [], structures_pdb: [], msa: [], pae: [], ptms: [], iptms: [], ipsaes: [], chainwise: [], chainwise_ipsaes: [], confidence: [], plddt: [], scores: [], raw: []]
    backends_out.each { backend, out ->
        combined.top_ranked_structure   << out.top_ranked_structure
        combined.structures             << out.structures
        combined.structures_pdb         << out.structures_pdb
        combined.msa                    << out.msa
        combined.pae                    << out.pae
        combined.ptms                   << out.ptms
        combined.iptms                  << out.iptms
        combined.ipsaes                 << out.ipsaes
        combined.chainwise             << out.chainwise
        combined.chainwise_ipsaes      << out.chainwise_ipsaes
        combined.confidence             << out.confidence
        combined.plddt                  << out.plddt
        combined.scores                 << out.scores
        combined.raw                    << out.raw
    }
    def ch_multiqc_metrics = Channel.empty()
    backends_out.each { backend, out ->
        ch_multiqc_metrics = ch_multiqc_metrics.mix(
            collectMultiqcMetrics([
                [ 'plddt', out.plddt ],
                [ 'msa',   out.msa ],
                [ 'ptm',   out.ptms ],
                [ 'iptm',  out.iptms ]
            ])
        )
    }

    // Mix the per-backend channels; inject instead of the spread operator,
    // which the Nextflow linter cannot parse. Optional emits can be empty
    // per backend; inject on an empty list yields Channel.empty().
    def ch_batch = [
        top_ranked_structure: combined.top_ranked_structure.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        structures          : combined.structures.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        structures_pdb      : combined.structures_pdb.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        msa                 : combined.msa.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        pae                 : combined.pae.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        ptms                : combined.ptms.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        iptms               : combined.iptms.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        ipsaes              : combined.ipsaes.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        chainwise           : combined.chainwise.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        chainwise_ipsaes    : combined.chainwise_ipsaes.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        confidence          : combined.confidence.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        plddt               : combined.plddt.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        scores              : combined.scores.inject(Channel.empty()) { acc, ch -> acc.mix(ch) },
        raw                 : combined.raw.inject(Channel.empty()) { acc, ch -> acc.mix(ch) }
    ]
    ch_top_ranked_structure = ch_batch.top_ranked_structure
    ch_structures = ch_batch.structures_pdb.mix(ch_batch.structures)

    emit:
    top_ranked_pdb  = ch_top_ranked_structure
    pdb             = ch_structures
    msa             = ch_batch.msa
    pae             = ch_batch.pae
    ptm             = ch_batch.ptms
    iptm            = ch_batch.iptms
    ipsae           = ch_batch.ipsaes
    chainwise       = ch_batch.chainwise
    chainwise_ipsae = ch_batch.chainwise_ipsaes
    confidence      = ch_batch.confidence
    plddt           = ch_batch.plddt
    multiqc_metrics = ch_multiqc_metrics
}
