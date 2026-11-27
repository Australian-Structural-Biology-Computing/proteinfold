process COLABFOLD2_INFERENCE {
    tag "$meta.id ($model_type)"
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
            if (dir.contains(':') || dir.contains(',')) {
                error("--colabfold2_compile_cache must not contain ':' or ',' (container bind delimiters)")
            }
            def quoted_dir = "'" + dir.replace("'", "'\\''") + "'"
            if (workflow.containerEngine in ['singularity', 'apptainer']) {
                "--bind ${quoted_dir}"
            } else if (workflow.containerEngine == 'docker') {
                "-v ${quoted_dir}:${quoted_dir}"
            }
        }
    }

    input:
    tuple val(meta), path(input), val(model_type), val(use_esm), path(params_path)
    val num_recycles

    output:
    path "raw/**"                                             , optional: true, emit: raw
    tuple val(meta), path("${meta.id}_colabfold2.*")          , optional: true, emit: top_ranked_structure
    tuple val(meta), path("raw/*_rank_*.cif")                 , optional: true, emit: structures
    tuple val(meta), path("raw/*_rank_*.pdb")                 , optional: true, emit: structures_pdb
    tuple val(meta), path("${meta.id}_colabfold2_msa.tsv")    , optional: true, emit: msa
    tuple val(meta), path("${meta.id}_1_pae.tsv")             , optional: true, emit: pae
    tuple val(meta), path("${meta.id}_ptm.tsv")               , optional: true, emit: ptms
    tuple val(meta), path("${meta.id}_iptm.tsv")              , optional: true, emit: iptms
    tuple val(meta), path("${meta.id}_ipsae.tsv")             , optional: true, emit: ipsaes
    tuple val(meta), path("${meta.id}_chainwise.tsv")         , optional: true, emit: chainwise
    tuple val(meta), path("${meta.id}_chainwise_ipsae.tsv")   , optional: true, emit: chainwise_ipsaes
    tuple val(meta), path("${meta.id}_plddt.tsv")             , optional: true, emit: plddt
    tuple val(meta), path("${meta.id}_confidence.tsv")        , optional: true, emit: confidence
    tuple val(meta), path("raw/*_scores_rank_*.json")         , optional: true, emit: scores
    tuple val("${task.process}"), val('colabfold'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"colabfold\"))' 2>/dev/null || echo \"unknown\""), emit: versions_colabfold, topic: versions
    tuple val("${task.process}"), val('colabfold_commit'), val('38683c21700408b6a29cf395e8ef07cf109fcb91'), emit: versions_colabfold_commit, topic: versions
    tuple val("${task.process}"), val('alphafold3_colabfold'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"alphafold3-colabfold\"))' 2>/dev/null || echo \"unknown\""), emit: versions_alphafold3_colabfold, topic: versions
    tuple val("${task.process}"), val('jax'), eval("python3 -c 'import jax; print(jax.__version__)' 2>/dev/null || echo \"unknown\""), emit: versions_jax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local COLABFOLD2_INFERENCE module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def use_esm_arg = use_esm ? '--use-esm' : ''
    def precision = task.ext.weights_precision ?: 'fp32'
    // Forward the global --random_seed; colabfold_batch defaults to 0 anyway.
    def random_seed_arg = params.random_seed != null ? "--random-seed ${params.random_seed}" : ''
    // Persistent JIT cache for the AF3-family backends: their pinned build has no
    // native compile-cache consumption, so jax's own persistent compilation cache
    // (JAX_COMPILATION_CACHE_DIR) covers them. The AlphaFold2 backends get
    // colabfold's per-shape executable cache via --compile-cache.
    // Same canonicalisation as the containerOptions bind above: the script
    // text carries the absolute path so the mount and the flags always agree.
    def compile_cache = params.colabfold2_compile_cache ? file(params.colabfold2_compile_cache).toAbsolutePath().toString() : ''
    def input_files = input instanceof Collection ? input.collect { file_path -> file_path.name } : [input.name]
    if (input_files.size() != 1) {
        error("COLABFOLD2_INFERENCE received multiple inputs for sample '${meta.id}'. Use --split_fasta to run each FASTA record independently.")
    }
    def require_msa_output = input_files.any { name -> name ==~ /(?i).*\.(fa|faa|fasta)$/ }
    // The alphafold2 backends write PDB structures; the AF3-family writes mmCIF.
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2')) ? 'pdb' : 'cif'
    def shell_quote = { name -> "'" + name.replace("'", "'\\''") + "'" }
    def quoted_input_files = input_files.collect(shell_quote).join(' ')
    def quoted_compile_cache = shell_quote.call(compile_cache)
    def compile_cache_arg = compile_cache ? "--compile-cache ${quoted_compile_cache}" : ''

    """
    if ! command -v colabfold_batch >/dev/null 2>&1; then
        echo "COLABFOLD2_INFERENCE: colabfold_batch is missing or not executable on the container PATH" >&2
        exit 1
    fi
    if ! python3 -c 'import colabfold, colabfold.alphafold3.predict, jax'; then
        echo "COLABFOLD2_INFERENCE: runtime import validation failed (colabfold/jax not importable, see traceback above)" >&2
        exit 1
    fi
    mkdir -p msa_in
    for f in ${quoted_input_files}; do
        ln -s "\$PWD/\$f" msa_in/
        if [[ "\$f" == *.json ]]; then
            python3 -c 'import json, sys; d=json.load(open(sys.argv[1])); n=len(d) if isinstance(d, list) else 1; sys.exit(f"ColabFold2 requires one AlphaFold3 job per samplesheet row; {sys.argv[1]} contains {n} jobs") if n != 1 else None' "\$f"
        fi
    done
    # alphafold2 backends consume ColabFold a3m queries; the AF3-family
    # consumes FASTA or AlphaFold3 JSON with embedded MSAs.

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
    export PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}"
    if [ ${quoted_compile_cache} != '' ]; then
        mkdir -p ${quoted_compile_cache}
        export JAX_COMPILATION_CACHE_DIR=${quoted_compile_cache}
    fi

    mkdir raw
    colabfold_batch_af3_chainwise.py \
        --model-type ${model_type} \
        ${use_esm_arg} \
        --calc-extra-ptm \
        --num-recycle ${num_recycles} \
        --weights-precision ${precision} \
        --data "${params_path}" \
        ${random_seed_arg} \
        ${compile_cache_arg} \
        $args \
        msa_in/ \
        raw/

    mapfile -t ranked_structures < <(find raw -maxdepth 1 -type f -name "*_rank_*.${structure_ext}" | sort)
    if (( \${#ranked_structures[@]} == 0 )); then
        echo "COLABFOLD2_INFERENCE produced no ranked structures" >&2
        exit 1
    fi
    cp "\${ranked_structures[0]}" "${meta.id}_colabfold2.${structure_ext}"

    mapfile -t score_files < <(find raw -maxdepth 1 -type f -name '*_scores_rank_*.json' | sort)
    if (( \${#score_files[@]} == 0 )); then
        echo "COLABFOLD2_INFERENCE produced no confidence JSON files" >&2
        exit 1
    fi

    extract_metrics.py --name "${meta.id}" \
        --colabfold_metrics_files "\${score_files[@]}" \
        --structs "\${ranked_structures[@]}"

    msa_file=\$(find raw -maxdepth 1 -type f -name '*.a3m' -print -quit)
    if [[ -n "\$msa_file" ]]; then
        extract_metrics.py --name "${meta.id}" --paired_a3m "\$msa_file"
        mv "${meta.id}_msa.tsv" "${meta.id}_colabfold2_msa.tsv"
    elif [[ "${require_msa_output}" == "true" ]]; then
        echo "COLABFOLD2_INFERENCE received FASTA input but the MSA server produced no a3m files" >&2
        exit 1
    else
        input_file=${quoted_input_files}
        if [[ "\$input_file" == *.a3m ]]; then
            extract_metrics.py --name "${meta.id}" --paired_a3m "\$input_file"
            mv "${meta.id}_msa.tsv" "${meta.id}_colabfold2_msa.tsv"
        elif [[ "\$input_file" == *.json ]]; then
            if af3_msa_json_to_a3m.py "\$input_file" --output "${meta.id}_input.a3m"; then
                extract_metrics.py --name "${meta.id}" --paired_a3m "${meta.id}_input.a3m"
                mv "${meta.id}_msa.tsv" "${meta.id}_colabfold2_msa.tsv"
            else
                conversion_status=\$?
                if [[ "\$conversion_status" -eq 2 || "\$conversion_status" -eq 3 ]]; then
                    touch "${meta.id}_colabfold2_msa.tsv"
                else
                    echo "COLABFOLD2_INFERENCE could not convert the input MSA" >&2
                    exit "\$conversion_status"
                fi
            fi
        else
            echo "COLABFOLD2_INFERENCE cannot extract an MSA from \$input_file" >&2
            exit 1
        fi
    fi

    touch \
        "${meta.id}_1_pae.tsv" \
        "${meta.id}_ptm.tsv" \
        "${meta.id}_iptm.tsv" \
        "${meta.id}_ipsae.tsv" \
        "${meta.id}_chainwise.tsv" \
        "${meta.id}_chainwise_ipsae.tsv" \
        "${meta.id}_plddt.tsv"
    cp "${meta.id}_plddt.tsv" "${meta.id}_confidence.tsv"

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
    """
    mkdir raw
    touch "${meta.id}_colabfold2.${structure_ext}"
    touch "raw/${meta.id}_rank_001_${model_type}_seed_000_sample_0.${structure_ext}"
    touch "raw/${meta.id}_scores_rank_001_${model_type}_seed_000_sample_0.json"
    touch "${meta.id}_colabfold2_msa.tsv"
    touch "${meta.id}_1_pae.tsv"
    touch "${meta.id}_ptm.tsv"
    touch "${meta.id}_iptm.tsv"
    touch "${meta.id}_ipsae.tsv"
    touch "${meta.id}_chainwise.tsv"
    touch "${meta.id}_chainwise_ipsae.tsv"
    touch "${meta.id}_plddt.tsv"
    touch "${meta.id}_confidence.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: 1.6.3
        colabfold_commit: 38683c21700408b6a29cf395e8ef07cf109fcb91
        alphafold3_colabfold: 3.1.18
        jax: unknown
    END_VERSIONS
    """
}
