include { batchTag } from '../../../subworkflows/local/utils_nfcore_proteinfold_pipeline'

def shell_quote(name) {
    "'" + name.replace("'", "'\\''") + "'"
}


process COLABFOLD2_INFERENCE_BATCH {
    tag "${batchTag(batch_meta)} ($model_type)"
    label 'process_medium'
    label 'process_gpu'

    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    // The compiled-kernel cache lives outside the work directory and is only
    // referenced inside the script, so container autoMounts never sees it:
    // bind it explicitly so users don't have to pass --singularity-args.
    // A bind source must exist before the container starts, which the
    // pipeline guarantees by creating it on the host at launch.
    containerOptions {
        if (params.colabfold2_compile_cache) {
            // Resolve against the launch directory and canonicalise: a bind
            // destination must be absolute (a relative --colabfold2_compile_cache
            // makes apptainer abort the container with 'destination must be an
            // absolute path').
            def dir = file(params.colabfold2_compile_cache).toAbsolutePath().toString()
            if (workflow.containerEngine in ['singularity', 'apptainer']) {
                "--bind ${dir}"
            } else if (workflow.containerEngine == 'docker') {
                "-v ${dir}:${dir}"
            }
        }
    }

    input:
    tuple val(batch_meta), path(query_files), val(model_type), val(use_esm), path(params_path)
    val num_recycles

    output:
    tuple val(batch_meta), path("*_colabfold2.*")              , emit: top_ranked_structures
    tuple val(batch_meta), path("raw/*_rank_*.cif")            , optional: true, emit: structures_cif
    tuple val(batch_meta), path("raw/*_rank_*.pdb")            , optional: true, emit: structures_pdb
    tuple val(batch_meta), path("*_colabfold2_metrics.jsonl.gz"), path("*_colabfold2_metrics.npz"), optional: true, emit: metrics_packed
    tuple val(batch_meta), path("raw/*.a3m")                   , optional: true, emit: msa_a3m
    path "raw/**"                                              , optional: true, emit: raw
    tuple val("${task.process}"), val('colabfold'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"colabfold\"))' 2>/dev/null || echo \"unknown\""), emit: versions_colabfold, topic: versions
    tuple val("${task.process}"), val('colabfold_commit'), val('38683c21700408b6a29cf395e8ef07cf109fcb91'), emit: versions_colabfold_commit, topic: versions
    tuple val("${task.process}"), val('alphafold3_colabfold'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"alphafold3-colabfold\"))' 2>/dev/null || echo \"unknown\""), emit: versions_alphafold3_colabfold, topic: versions
    tuple val("${task.process}"), val('jax'), eval("python3 -c 'import jax; print(jax.__version__)' 2>/dev/null || echo \"unknown\""), emit: versions_jax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local COLABFOLD2_INFERENCE_BATCH module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def use_esm_arg = use_esm ? '--use-esm' : ''
    def precision = task.ext.weights_precision ?: 'fp32'
    // Forward the global --random_seed; colabfold_batch defaults to 0 anyway.
    def random_seed_arg = params.random_seed != null ? "--random-seed ${params.random_seed}" : ''
    // Persistent JIT cache for the AF3-family backends: their pinned build has no
    // native compile-cache consumption, so jax's own persistent compilation cache
    // (JAX_COMPILATION_CACHE_DIR) covers them. The AlphaFold2 backends get
    // colabfold's per-shape executable cache via the --compile-cache ext.args entry.
    // Same canonicalisation as the containerOptions bind above: the script
    // text carries the absolute path so the mount and the flags always agree.
    def compile_cache = params.colabfold2_compile_cache ? file(params.colabfold2_compile_cache).toAbsolutePath().toString() : ''
    // The alphafold2 backends write PDB structures; the AF3-family writes mmCIF.
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2')) ? 'pdb' : 'cif'
    def files = query_files instanceof Collection ? query_files.collect() : [query_files]
    def metas = batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]
    if (metas.size() != files.size()) {
        error("COLABFOLD2_INFERENCE_BATCH received ${metas.size()} sample metadata entries but ${files.size()} input files")
    }
    // colabfold_batch sanitises job names to [A-Za-z0-9_.-]; stage every
    // query under the sanitised sample id so outputs stay addressable.
    def sample_ids = metas.collect { (it.id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    if (sample_ids.toSet().size() != sample_ids.size()) {
        error("COLABFOLD2_INFERENCE_BATCH received sample IDs that collide after colabfold job-name sanitisation: ${sample_ids}")
    }
    sample_ids.each { id ->
        if (id =~ /[\t\n\r]/) {
            error("Sample ID '${id}' contains characters that cannot be staged for batched inference")
        }
    }
    def stage_cmds = [metas, files].transpose().collect { meta, file ->
        def name = file.name
        def dot = name.lastIndexOf('.')
        def ext = dot >= 0 ? name.substring(dot) : ''
        "cp -- ${shell_quote(name)} ${shell_quote("msa_in/${(meta.id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')}${ext}")}"
    }
    def ids_content = sample_ids.collect { shell_quote(it) }.join(' ')

    // Per-sample metrics are packed into two per-batch stores (see
    // bin/pack_metrics.py); per-sample metric TSVs are never produced.
    def per_id_cmds = """
    [ -f "raw/\$id.a3m" ] || touch "raw/\$id.a3m"
"""
    def post_loop_cmds = """

pack_metrics.py --ids sample_ids.txt --struct-ext ${structure_ext} --raw-dir raw \
    --out-scalars '${sample_ids[0]}_colabfold2_metrics.jsonl.gz' \
    --out-arrays '${sample_ids[0]}_colabfold2_metrics.npz'
"""
    // The packed stores carry everything downstream consumers need. Delete
    // colabfold's bulky per-rank scores JSONs and the msa_in staging copies
    // once packed so they do not weigh down the work directory -- unless
    // --save_intermediates was requested, which publishes raw/ verbatim.
    def raw_cleanup_cmds = params.save_intermediates ? '' : """
    rm -f raw/*_scores_rank_*.json
    rm -rf msa_in __pycache__
"""
    """
if ! command -v colabfold_batch >/dev/null 2>&1; then
    echo "COLABFOLD2_INFERENCE_BATCH: colabfold_batch is missing or not executable on the container PATH" >&2
    exit 1
fi
if ! python3 -c 'import colabfold, colabfold.alphafold3.predict, jax'; then
    echo "COLABFOLD2_INFERENCE_BATCH: runtime import validation failed (colabfold/jax not importable, see traceback above)" >&2
    exit 1
fi
mkdir -p msa_in
${stage_cmds.join('\n')}
# Write the sanitised sample ids one per line (printf keeps this robust
# against the template's uneven indentation, unlike a shell heredoc).
printf '%s\\n' ${ids_content} > sample_ids.txt
python3 - <<'PY'
import glob, json, os, sys
for path in sorted(glob.glob('msa_in/*.json')):
    with open(path) as handle:
        document = json.load(handle)
    if isinstance(document, list):
        if len(document) != 1:
            sys.exit(f"{path} contains {len(document)} AlphaFold3 jobs; ColabFold2 requires one job per samplesheet row")
        job = document[0]
    else:
        job = document
    # Pin the colabfold job name to the sample id so the outputs of a
    # batched run stay addressable per sample.
    job['name'] = os.path.splitext(os.path.basename(path))[0]
    with open(path, 'w') as handle:
        json.dump(document, handle, indent=2)
PY
for f in msa_in/*.fa msa_in/*.faa msa_in/*.fasta; do
    [ -e "\$f" ] || continue
    nrecords=\$(grep -c '^>' "\$f" || true)
    if [ "\$nrecords" != "1" ]; then
        echo "COLABFOLD2_INFERENCE_BATCH: \$f contains \$nrecords FASTA records; use --split_fasta to run each record independently" >&2
        exit 1
    fi
done

cat > sitecustomize.py <<'EOF'
def _register_missing_gpu_kinds():
    try:
        from jax._src.pallas.triton import gpu_info as gi
    except Exception:
        return
    def add(kind, arch, cc):
        if kind not in gi.registry:
            gi.registry[kind] = lambda: gi.GpuInfo(gpu_version=None, arch_name=arch, compute_capability=cc)
    add("NVIDIA L40S", "8.9", 89)
    add("NVIDIA A10G", "8.6", 86)
    for kind in ("Tesla V100-SXM2-32GB", "Tesla V100-SXM2-16GB",
                 "Tesla V100-PCIE-32GB", "Tesla V100-PCIE-16GB"):
        add(kind, "7.0", 70)
_register_missing_gpu_kinds()
EOF
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}"
if [ '${compile_cache}' != '' ]; then
    mkdir -p '${compile_cache}'
    export JAX_COMPILATION_CACHE_DIR='${compile_cache}'
fi

mkdir raw
colabfold_batch_af3_chainwise.py \\
    --model-type ${model_type} \\
    ${use_esm_arg} \\
    --calc-extra-ptm \\
    --num-recycle ${num_recycles} \\
    --weights-precision ${precision} \\
    --data "${params_path}" \\
    ${random_seed_arg} \\
    $args \\
    msa_in/ \\
    raw/

while IFS= read -r id; do
    mapfile -t ranked_structures < <(find raw -maxdepth 1 -type f -name "\${id}_rank_*.${structure_ext}" | sort)
    if (( \${#ranked_structures[@]} == 0 )); then
        echo "COLABFOLD2_INFERENCE_BATCH produced no ranked structures for \${id}" >&2
        exit 1
    fi
    cp "\${ranked_structures[0]}" "\${id}_colabfold2.${structure_ext}"

    mapfile -t score_files < <(find raw -maxdepth 1 -type f -name "\${id}_scores_rank_*.json" | sort)
    if (( \${#score_files[@]} == 0 )); then
        echo "COLABFOLD2_INFERENCE_BATCH produced no confidence JSON files for \${id}" >&2
        exit 1
    fi
${per_id_cmds}
done < sample_ids.txt
${post_loop_cmds}
${raw_cleanup_cmds}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    colabfold: \$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("colabfold"))' 2>/dev/null || echo unknown)
    colabfold_commit: 38683c21700408b6a29cf395e8ef07cf109fcb91
    alphafold3_colabfold: \$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("alphafold3-colabfold"))' 2>/dev/null || echo unknown)
    jax: \$(python3 -c 'import jax; print(jax.__version__)' 2>/dev/null || echo unknown)
END_VERSIONS
    """

    stub:
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2')) ? 'pdb' : 'cif'
    def stub_metas = batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]
    def first_id = (stub_metas[0].id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
    def stub_commands = stub_metas.collect { meta -> """
    touch "${meta.id}_colabfold2.${structure_ext}"
    touch "raw/${meta.id}_rank_001_${model_type}_seed_000_sample_0.${structure_ext}"
    touch "raw/${meta.id}_scores_rank_001_${model_type}_seed_000_sample_0.json"
    touch "raw/${meta.id}.a3m"
    """ }.join('\n')
    """
    mkdir raw
${stub_commands}
    touch "${first_id}_colabfold2_metrics.jsonl.gz" "${first_id}_colabfold2_metrics.npz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: 1.6.3
        colabfold_commit: 38683c21700408b6a29cf395e8ef07cf109fcb91
        alphafold3_colabfold: 3.1.18
        jax: unknown
    END_VERSIONS
    """
}
