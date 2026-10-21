process MODELCIF_ASSEMBLE {
    tag "$meta.id-$meta.model"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/keiran-rowell-unsw/proteinfold-modelcif:1.7'

    // Placeholder TSVs keep the tuple complete; populate_modelcif.py omits them from metadata.
    // TODO: Unsure at this stage best way forward with chain-wise paths. At least modelCIF has the asymmetric unit entity, so we *could* DUMMY_FILE, but I do like the philosophy of don't check things that *can't* exist
    input:
    tuple val(meta), path(structs), path(msa), path(plddt), path(pae), path(ptm), path(iptm), path(versions_yml, stageAs: 'versions_input/versions.yml'), path(software_details)

    output:
    tuple val(meta), path("*.{mmcif,bcif}"), emit: modelcif
    path "pae/*.tsv", optional: true, emit: pae
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    // This will take a maximalist approach to populating the modelCIF
    // Every single structure from a sequence prediction method will be captured in .ModelGroup
    // Every common metric will be captured as a modelCIF qa_metric
    // The advantage is all software metadata and protocol steps are shared between all structures in the .ModelGroup
    // BinaryCIF output: pass --write_binary via task.ext.args.
    script:
    def args = task.ext.args ?: ''
    def no_template_search = meta.no_template_search ? '--no-template-search' : ''
    def missing_metrics = meta.missing_metrics ? "--missing-metrics ${meta.missing_metrics.join(',')}" : ''
    def prediction_container = meta.prediction_container ? "--container_image ${meta.prediction_container}" : ''
    def stage_pae = meta.missing_metrics?.contains('pae') ? '' : "mkdir -p pae && cp ${pae} pae/"
    def pae_inputs = meta.missing_metrics?.contains('pae') ? "${pae}" : 'pae/*.tsv'
    def modelcif_params = (meta.modelcif_params ?: [:]).collect { key, value -> "--param ${key}=${value}" }.join(' ')
    def plddt_scale = meta.plddt_scale ? "--plddt-scale ${meta.plddt_scale}" : ''
    def output_prefix = "${meta.id}_${meta.model}"
    """
    ${stage_pae}

    populate_modelcif.py \\
        --structs ${structs} \\
        --msa ${msa} \\
        --plddt ${plddt} \\
        --pae ${pae_inputs} \\
        --ptm ${ptm} \\
        --iptm ${iptm} \\
        --name ${meta.id} \\
        --prog ${meta.model} \\
        --versions_yml ${versions_yml} \\
        --msa_tool ${meta.msa_tool ?: 'None'} \\
        --software_details ${software_details} \\
        ${prediction_container} \\
        ${missing_metrics} \\
        ${modelcif_params} \\
        ${plddt_scale} \\
        ${no_template_search} \\
        $args

    python3 -c 'import glob, modelcif.reader; outputs = [f for f in glob.glob("${output_prefix}.mmcif") + glob.glob("${output_prefix}.bcif")]; assert len(outputs) == 1, f"Expected one ModelCIF output, found: {outputs}"; output = outputs[0]; handle = open(output, "rb" if output.endswith(".bcif") else "r"); systems = modelcif.reader.read(handle); handle.close(); assert systems, f"No ModelCIF systems found in {output}"'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        modelcif: \$(python3 -c "import modelcif; print(modelcif.__version__)" 2>/dev/null || echo "unknown")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

}
