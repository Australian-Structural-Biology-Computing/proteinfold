process COLABFOLD2_INFERENCE {
    tag "$meta.id ($model_type)"
    label 'process_medium'
    label 'process_gpu'

    container "ghcr.io/jscgh/proteinfold_colabfold:af3-preview-e9e8bd0"

    input:
    tuple val(meta), path(input)
    val model_type
    val use_esm
    val public_mode
    val num_recycles
    path params_path

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
    tuple val(meta), path("${meta.id}_chainwise_iptm.tsv")    , optional: true, emit: chainwise_iptms
    tuple val(meta), path("${meta.id}_chainwise_ipsae.tsv")   , optional: true, emit: chainwise_ipsaes
    tuple val(meta), path("${meta.id}_plddt.tsv")             , optional: true, emit: plddt
    tuple val(meta), path("${meta.id}_confidence.tsv")        , optional: true, emit: confidence
    tuple val(meta), path("raw/*_scores_rank_*.json")         , optional: true, emit: scores
    tuple val("${task.process}"), val('colabfold'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"colabfold\"))' 2>/dev/null || echo \"unknown\""), emit: versions_colabfold, topic: versions
    tuple val("${task.process}"), val('colabfold_commit'), val('e9e8bd0b8b6c26f9865808070342d0d785dada2e'), emit: versions_colabfold_commit, topic: versions
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
    def input_files = input instanceof Collection ? input.collect { file_path -> file_path.name } : [input.name]
    if (input_files.size() != 1) {
        error("COLABFOLD2_INFERENCE received multiple inputs for sample '${meta.id}'. Use --split_fasta to run each FASTA record independently.")
    }
    def require_msa_output = input_files.any { name -> name ==~ /(?i).*\.(fa|faa|fasta)$/ }
    // The alphafold2 backends write PDB structures; the AF3-family writes mmCIF.
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2') || model_type.startsWith('deepfold')) ? 'pdb' : 'cif'
    def shell_quote = { name -> "'" + name.replace("'", "'\\''") + "'" }
    def quoted_input_files = input_files.collect(shell_quote).join(' ')

    """
    if ! command -v colabfold_batch >/dev/null 2>&1; then
        echo "COLABFOLD2_INFERENCE: colabfold_batch is missing or not executable on the container PATH" >&2
        exit 1
    fi
    if ! python3 -c 'import colabfold, jax'; then
        echo "COLABFOLD2_INFERENCE: runtime import validation failed (colabfold/jax not importable, see traceback above)" >&2
        exit 1
    fi
    mkdir -p msa_in
    for f in ${quoted_input_files}; do
        if [[ "\$f" == *.json ]]; then
            # copy so the job-name normalisation below cannot touch the staged input
            cp -- "\$PWD/\$f" msa_in/
        else
            ln -s "\$PWD/\$f" msa_in/
        fi
    done
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
        # Pin the colabfold job name to the sample id so the outputs are
        # addressable by sample regardless of the name inside the JSON.
        job['name'] = os.path.splitext(os.path.basename(path))[0]
        with open(path, 'w') as handle:
            json.dump(document, handle, indent=2)
    PY
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

    mkdir raw
    colabfold_batch \
        --model-type ${model_type} \
        ${use_esm_arg} \
        --num-recycle ${num_recycles} \
        --weights-precision ${precision} \
        --data "${params_path}" \
        ${random_seed_arg} \
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
                if [[ "\$conversion_status" -eq 2 ]]; then
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
        "${meta.id}_chainwise_iptm.tsv" \
        "${meta.id}_chainwise_ipsae.tsv" \
        "${meta.id}_plddt.tsv"
    cp "${meta.id}_plddt.tsv" "${meta.id}_confidence.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: \$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("colabfold"))' 2>/dev/null || echo unknown)
        colabfold_commit: e9e8bd0b8b6c26f9865808070342d0d785dada2e
        alphafold3_colabfold: \$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("alphafold3-colabfold"))' 2>/dev/null || echo unknown)
        jax: \$(python3 -c 'import jax; print(jax.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2') || model_type.startsWith('deepfold')) ? 'pdb' : 'cif'
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
    touch "${meta.id}_chainwise_iptm.tsv"
    touch "${meta.id}_chainwise_ipsae.tsv"
    touch "${meta.id}_plddt.tsv"
    touch "${meta.id}_confidence.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: 1.6.3
        colabfold_commit: e9e8bd0b8b6c26f9865808070342d0d785dada2e
        alphafold3_colabfold: 3.1.18
        jax: unknown
    END_VERSIONS
    """
}

process COLABFOLD2_INFERENCE_BATCH {
    tag "batch_${batch_meta*.id.join('_')} ($model_type)"
    label 'process_medium'
    label 'process_gpu'

    container "ghcr.io/jscgh/proteinfold_colabfold:af3-preview-e9e8bd0"

    input:
    tuple val(batch_meta), path(query_files)
    val model_type
    val use_esm
    val public_mode
    val num_recycles
    path params_path

    output:
    tuple val(batch_meta), path("*_colabfold2.*")              , emit: top_ranked_structures
    tuple val(batch_meta), path("raw/*_rank_*.cif")            , optional: true, emit: structures_cif
    tuple val(batch_meta), path("raw/*_rank_*.pdb")            , optional: true, emit: structures_pdb
    tuple val(batch_meta), path("*_colabfold2_msa.tsv")        , emit: msa
    tuple val(batch_meta), path("*_1_pae.tsv")                 , optional: true, emit: pae
    tuple val(batch_meta), path("*_ptm.tsv")                   , optional: true, emit: ptms
    tuple val(batch_meta), path("*_iptm.tsv")                  , optional: true, emit: iptms
    tuple val(batch_meta), path("*_ipsae.tsv")                 , optional: true, emit: ipsaes
    tuple val(batch_meta), path("*_chainwise_iptm.tsv")        , optional: true, emit: chainwise_iptms
    tuple val(batch_meta), path("*_chainwise_ipsae.tsv")       , optional: true, emit: chainwise_ipsaes
    tuple val(batch_meta), path("*_plddt.tsv")                 , emit: plddt
    tuple val(batch_meta), path("*_confidence.tsv")            , emit: confidence
    tuple val(batch_meta), path("raw/*_scores_rank_*.json")    , optional: true, emit: scores
    path "raw/**"                                              , optional: true, emit: raw
    tuple val("${task.process}"), val('colabfold'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"colabfold\"))' 2>/dev/null || echo \"unknown\""), emit: versions_colabfold, topic: versions

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
    // The alphafold2 backends write PDB structures; the AF3-family writes mmCIF.
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2') || model_type.startsWith('deepfold')) ? 'pdb' : 'cif'
    def shell_quote = { name -> "'" + name.replace("'", "'\\''") + "'" }
    def files = query_files instanceof Collection ? query_files.collect() : [query_files]
    def metas = batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]
    if (metas.size() != files.size()) {
        error("COLABFOLD2_INFERENCE_BATCH received ${metas.size()} sample metadata entries but ${files.size()} input files")
    }
    def sample_ids = metas.collect { it.id as String }
    if (sample_ids.toSet().size() != sample_ids.size()) {
        error("COLABFOLD2_INFERENCE_BATCH received duplicate sample IDs: ${sample_ids}")
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
        "cp -- ${shell_quote(name)} ${shell_quote("msa_in/${meta.id}${ext}")}"
    }
    def ids_content = sample_ids.join('\n')

    """
    if ! command -v colabfold_batch >/dev/null 2>&1; then
        echo "COLABFOLD2_INFERENCE_BATCH: colabfold_batch is missing or not executable on the container PATH" >&2
        exit 1
    fi
    if ! python3 -c 'import colabfold, jax'; then
        echo "COLABFOLD2_INFERENCE_BATCH: runtime import validation failed (colabfold/jax not importable, see traceback above)" >&2
        exit 1
    fi
    mkdir -p msa_in
    ${stage_cmds.join('\n')}
    cat > sample_ids.txt <<'COLABFOLD2_BATCH_IDS'
    ${ids_content}
    COLABFOLD2_BATCH_IDS
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
    export PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}"

    mkdir raw
    colabfold_batch \\
        --model-type ${model_type} \\
        ${use_esm_arg} \\
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

        extract_metrics.py --name "\$id" \\
            --colabfold_metrics_files "\${score_files[@]}" \\
            --structs "\${ranked_structures[@]}"

        msa_file=\$(find raw -maxdepth 1 -type f -name "\${id}.a3m" -print -quit)
        if [[ -n "\$msa_file" ]]; then
            extract_metrics.py --name "\$id" --paired_a3m "\$msa_file"
            mv "\${id}_msa.tsv" "\${id}_colabfold2_msa.tsv"
        else
            touch "\${id}_colabfold2_msa.tsv"
        fi

        for metric in 1_pae ptm iptm ipsae chainwise_iptm chainwise_ipsae plddt; do
            [ -f "\${id}_\${metric}.tsv" ] || touch "\${id}_\${metric}.tsv"
        done
        cp "\${id}_plddt.tsv" "\${id}_confidence.tsv"
    done < sample_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: \$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("colabfold"))' 2>/dev/null || echo unknown)
        colabfold_commit: e9e8bd0b8b6c26f9865808070342d0d785dada2e
        alphafold3_colabfold: \$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("alphafold3-colabfold"))' 2>/dev/null || echo unknown)
        jax: \$(python3 -c 'import jax; print(jax.__version__)' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    def structure_ext = (model_type == 'auto' || model_type.startsWith('alphafold2') || model_type.startsWith('deepfold')) ? 'pdb' : 'cif'
    def stub_metas = batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]
    def stub_commands = stub_metas.collect { meta -> """
    touch "${meta.id}_colabfold2.${structure_ext}"
    touch "raw/${meta.id}_rank_001_${model_type}_seed_000_sample_0.${structure_ext}"
    touch "raw/${meta.id}_scores_rank_001_${model_type}_seed_000_sample_0.json"
    touch "${meta.id}_colabfold2_msa.tsv"
    touch "${meta.id}_1_pae.tsv"
    touch "${meta.id}_ptm.tsv"
    touch "${meta.id}_iptm.tsv"
    touch "${meta.id}_ipsae.tsv"
    touch "${meta.id}_chainwise_iptm.tsv"
    touch "${meta.id}_chainwise_ipsae.tsv"
    touch "${meta.id}_plddt.tsv"
    touch "${meta.id}_confidence.tsv"
    """ }.join('\n')
    """
    mkdir raw
${stub_commands}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold: 1.6.3
        colabfold_commit: e9e8bd0b8b6c26f9865808070342d0d785dada2e
        alphafold3_colabfold: 3.1.18
        jax: unknown
    END_VERSIONS
    """
}
