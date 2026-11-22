/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified proteinfold workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Single workflow for the pipeline: users choose one MSA generation module
    (--msa_tool: colabfold_search | alphafold2_msa | alphafold3_msa, or
    --use_msa_server for the ColabFold webserver), and every structure
    prediction step runs through the colabfold2 backends.
*/

include { COLABFOLD2_INFERENCE_BATCH as ALPHAFOLD2 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as ALPHAFOLD3 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as BOLTZ2 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as OPENFOLD3 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as OPENBIND0 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as PROTEINX2 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as CHAI1 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as INTELLIFOLD2 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as OPENDDE } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as ROSETTAFOLD3 } from '../modules/local/colabfold2_inference'
include { COLABFOLD2_INFERENCE_BATCH as ESMFOLD2 } from '../modules/local/colabfold2_inference'

include { TO_COLABFOLD_FASTA } from '../modules/local/to_colabfold_fasta'
include { ALPHAFOLD3_JSON_TO_COLABFOLD_A3M   } from '../modules/local/af3_msa_json_to_a3m'
include { MSA_COLABFOLD_SEARCH } from '../subworkflows/local/msa_colabfold_search'
include { MSA_ALPHAFOLD2       } from '../subworkflows/local/msa_alphafold2'
include { MSA_ALPHAFOLD3       } from '../subworkflows/local/msa_alphafold3'


// Mirrors colabfold's safe_filename: colabfold_batch sanitises query/job
// names to [A-Za-z0-9_.-], so output files are prefixed with the sanitised
// sample id rather than the raw one.
def safeColabfoldName(name) {
    (name as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

// Map a file produced by COLABFOLD2_INFERENCE_BATCH to its sample by matching
// the longest sanitised sample-id prefix (sample ids may be prefixes of one
// another).
def matchBatchSample(batch_metas, filename) {
    // The batch lists travel through nextflow operators that may share or
    // reuse their buffers across concurrent consumers; iterate over private
    // copies so a concurrent structural modification cannot break the pass.
    def ids = new ArrayList(batch_metas).collect { safeColabfoldName(it.id) }.sort { -it.length() }
    def id = ids.find { filename.startsWith(it + '_') || filename.startsWith(it + '.') }
    if (id == null) {
        error("COLABFOLD2_INFERENCE_BATCH produced file '${filename}' that does not match any sample in the batch (${ids.join(', ')})")
    }
    return new ArrayList(batch_metas).find { safeColabfoldName(it.id) == id }
}

// One file per sample: emit [meta, file] tuples.
def splitBatchSingleFiles(out_channel) {
    out_channel.flatMap { batch_meta, files ->
        def list = files instanceof List ? files : [files]
        list.collect { file -> [matchBatchSample(batch_meta, file.name), file] }
    }
}

// Batched a3m outputs: one [meta, file] tuple per sample.
def splitBatchA3M(out_channel) {
    out_channel.flatMap { batch_meta, files ->
        def list = files instanceof List ? files : [files]
        list.collect { file -> [matchBatchSample(batch_meta, file.name), file] }
    }
}

// Several files per sample (ranked structures): emit [meta, [files]] sorted by rank.
def splitBatchGroupedFiles(out_channel) {
    out_channel.flatMap { batch_meta, files ->
        def list = files instanceof List ? files : [files]
        def grouped = list
            .collect { file -> [matchBatchSample(batch_meta, file.name), file] }
            .groupBy { it[0] }
        // Private copy: the batch metadata list is shared across concurrent
        // consumers of the same task output.
        new ArrayList(batch_meta).collect { meta ->
            def sample_files = (grouped[meta] ?: []).collect { it[1] }.sort { file ->
                def matcher = file.name =~ /_rank_(\d+)_/
                matcher.find() ? matcher.group(1).toInteger() : Integer.MAX_VALUE
            }
            [meta, sample_files]
        }
    }
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
                tuple(new ArrayList(batch*.getAt(0)), new ArrayList(batch*.getAt(1)), model_types[0], use_esms[0], dirs[0])
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
            // JSON inputs are reduced to their protein chains first. One
            // conversion task converts every JSON sample and keeps one query
            // FASTA per sample for the per-sample inference.
            TO_COLABFOLD_FASTA(
                ch_input_by_ext.json
                    .collect(flat: false)
                    .filter { rows -> !rows.isEmpty() }
                    .map { rows -> [rows, rows.collect { it[1] }.unique { it.toString() }] },
                1,
                '',
                true
            )
            ch_af2_inference_input = ch_input_by_ext.fasta
                .mix(
                    TO_COLABFOLD_FASTA.out.queries
                        .flatMap { pairs, files ->
                            (files instanceof List ? files : [files]).collect { file ->
                                def sid = file.name.replaceAll(/\.fasta$/, '')
                                def pair = pairs.find { (it[0].id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_') == sid }
                                pair == null ? error("Converted query ${file.name} does not match any input sample") : [pair[0], file]
                            }
                        }
                )
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

    // Inference batching: N samples share one colabfold_batch task (one
    // model/weight load for the whole batch); 1 runs one sample per task.
    // Every backend runs the single batched process, for batch size 1 too.
    def inference_batch = inference_batch_size instanceof Number ?
        inference_batch_size : inference_batch_size.toString().toInteger()
    if (inference_batch == 0 || inference_batch < -1) {
        error("Inference batch size must be -1 or a positive integer, got ${inference_batch}")
    }
    def effective_inference_batch = msa_tool == 'msa_server' ? 1 : inference_batch
    if (msa_tool == 'msa_server' && inference_batch != 1) {
        log.warn "Inference batching is not supported with --use_msa_server; running inference per sample."
    }

    def ch_backend_inputs = ch_backends.collectEntries { backend, ch ->
        [backend, batchInferenceInputs(ch, effective_inference_batch)]
    }

    def backends_out = [:]
    if (af2_specs) {
        ALPHAFOLD2(ch_backend_inputs['ALPHAFOLD2'], num_recycles)
        backends_out['ALPHAFOLD2'] = ALPHAFOLD2.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'ALPHAFOLD3' }) {
        ALPHAFOLD3(ch_backend_inputs['ALPHAFOLD3'], num_recycles)
        backends_out['ALPHAFOLD3'] = ALPHAFOLD3.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'BOLTZ2' }) {
        BOLTZ2(ch_backend_inputs['BOLTZ2'], num_recycles)
        backends_out['BOLTZ2'] = BOLTZ2.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'OPENFOLD3' }) {
        OPENFOLD3(ch_backend_inputs['OPENFOLD3'], num_recycles)
        backends_out['OPENFOLD3'] = OPENFOLD3.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'OPENBIND0' }) {
        OPENBIND0(ch_backend_inputs['OPENBIND0'], num_recycles)
        backends_out['OPENBIND0'] = OPENBIND0.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'PROTEINX2' }) {
        PROTEINX2(ch_backend_inputs['PROTEINX2'], num_recycles)
        backends_out['PROTEINX2'] = PROTEINX2.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'CHAI1' }) {
        CHAI1(ch_backend_inputs['CHAI1'], num_recycles)
        backends_out['CHAI1'] = CHAI1.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'INTELLIFOLD2' }) {
        INTELLIFOLD2(ch_backend_inputs['INTELLIFOLD2'], num_recycles)
        backends_out['INTELLIFOLD2'] = INTELLIFOLD2.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'OPENDDE' }) {
        OPENDDE(ch_backend_inputs['OPENDDE'], num_recycles)
        backends_out['OPENDDE'] = OPENDDE.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'ROSETTAFOLD3' }) {
        ROSETTAFOLD3(ch_backend_inputs['ROSETTAFOLD3'], num_recycles)
        backends_out['ROSETTAFOLD3'] = ROSETTAFOLD3.out
    }
    if (mode_specs.any { COLABFOLD2_BACKENDS[it.model_type] == 'ESMFOLD2' }) {
        ESMFOLD2(ch_backend_inputs['ESMFOLD2'], num_recycles)
        backends_out['ESMFOLD2'] = ESMFOLD2.out
    }

    // Split every backend's batch outputs back into per-sample streams and
    // group the packed metric stores per public mode for MultiQC. The metas
    // already carry the public mode and MSA source (tagged at dispatch), so
    // publication stays collision-free without further tagging.
    def ch_top_ranked_structure = Channel.empty()
    def ch_structures_cif = Channel.empty()
    def ch_structures_pdb = Channel.empty()
    def ch_msa = Channel.empty()
    def ch_packed_metrics = Channel.empty()
    backends_out.each { backend, out ->
        ch_top_ranked_structure = ch_top_ranked_structure.mix(splitBatchSingleFiles(out.top_ranked_structures))
        ch_structures_cif       = ch_structures_cif.mix(splitBatchGroupedFiles(out.structures_cif))
        ch_structures_pdb       = ch_structures_pdb.mix(splitBatchGroupedFiles(out.structures_pdb))
        ch_msa                  = ch_msa.mix(splitBatchA3M(out.msa_a3m))
        ch_packed_metrics       = ch_packed_metrics.mix(out.metrics_packed)
    }
    def ch_structures = ch_structures_pdb.mix(ch_structures_cif)

    // One grouped task per public mode receives every batch's scalar + array
    // store; the MultiQC contents generator reads the stores natively.
    def ch_multiqc_metrics = ch_packed_metrics
        .map { batch_metas, scalars, arrays ->
            def mode = (batch_metas instanceof List ? batch_metas[0] : batch_metas).model
            [[id: mode, model: mode], [scalars, arrays]]
        }
        .groupTuple()
        .map { mode_meta, file_lists -> [mode_meta, file_lists.flatten()] }

    emit:
    top_ranked_pdb  = ch_top_ranked_structure
    pdb             = ch_structures
    msa             = ch_msa
    multiqc_metrics = ch_multiqc_metrics
    metrics_packed  = ch_packed_metrics
}
