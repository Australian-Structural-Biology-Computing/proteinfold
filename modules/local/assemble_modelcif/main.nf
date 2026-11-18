process MODELCIF_ASSEMBLE_BATCH {
    tag   "batch_${batch_meta*.id.join('_')}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/keiran-rowell-unsw/proteinfold-modelcif:1.7'

    input:
    tuple val(batch_meta), path(structures), path(msa_files), path(metrics_stores), val(manifest)
    path(versions_yml, stageAs: 'versions_input/versions.yml')
    path(software_details)

    output:
    tuple val(batch_meta), path("*.{mmcif,bcif}"), emit: modelcif
    path "pae/*.npz", optional: true, emit: pae
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Assembles ModelCIF for a whole batch: the driver reads the packed
    // metric stores natively and derives each sample's missing metrics from
    // the store contents. Template written at column 0 (see
    // colabfold2_inference for why).
    def args = task.ext.args ?: ''
    def batch_head = batch_meta instanceof Collection ? batch_meta[0] : batch_meta
    def model = batch_head.mode_short ?: batch_head.model
    // The alphafold2 backends write PDB structures; the AF3-family writes mmCIF.
    def struct_ext = (model == 'alphafold2') ? 'pdb' : 'cif'
    def msa_tool = batch_head.msa_tool ?: 'None'
    """
cat > sample_manifest.tsv <<BATCH_MANIFEST
${manifest}
BATCH_MANIFEST

export PYTHONDONTWRITEBYTECODE=1

populate_modelcif.py \\
    --manifest sample_manifest.tsv \\
    --struct-ext ${struct_ext} \\
    --prog ${model} \\
    --versions_yml versions_input/versions.yml \\
    --software_details ${software_details} \\
    --msa_tool ${msa_tool} ${args}

# ...then rewrite versions.yml in the collected-versions topic format.
cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        modelcif: \$(python3 -c "import modelcif; print(modelcif.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
"""

    stub:
    def stub_model = (batch_meta instanceof Collection ? batch_meta[0] : batch_meta).with { m -> m.mode_short ?: m.model } ?: 'colabfold2'
    def stub_files = (batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]).collect { meta -> """
    touch ${meta.id}_${meta.mode_short}_${meta.msa_short}.mmcif
    """ }.join('\n')
    """
${stub_files}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        modelcif: \$(python3 -c "import modelcif; print(modelcif.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
"""
}
