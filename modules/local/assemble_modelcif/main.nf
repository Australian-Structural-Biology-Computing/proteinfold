process ASSEMBLE_MODELCIF {
    tag "$meta.id-$meta.model"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "python:3.12-slim"

    // Unpack the tuple assuming every value used is at the path. Pass along a DUMMY_FILE so that the path is always occupied.
    // The populate_modelcif.py can explicitly trigger if the arg.metric is a DUMMY_FILE, and if so, appropriately account for the non-existence of that in the model metadata
    // This is also means no hand-crafted exceptions for particular programs that doesn't make use of some input (MSA - EMSFold), they all obey the same take: pipleline logic
    // TODO: Unsure at this stage best way forward with chain-wise paths. At least modelCIF has the asymmetric unit entity, so we *could* DUMMY_FILE, but I do like the philosophy of don't check things that *can't* exist
    input:
    tuple val(meta), path(structs), path(msa), path(plddt), path(pae), path(ptm), path(iptm), path(versions_yml)
    // TODO: --plddt-scale covers PLDDT/PLDDT01/PLDDTAllAtom/PLDDTAllAtom01; still need per-program scale selection
    // and removal of the averaging in EXTRACT_METRICS
    // TODO: A space will be made for path(ipsae) once 1) it's captured 2) an ipsae custom class extends the modelCIF construction
    // Rank mapping is guaranteed: extract_metrics.py emits rank_0..N; populate_modelcif.py rejects TSVs without them.
    // Seed: --seed takes a {name}_seed.tsv (see assets/DUMMY_SEED.tsv) or an int; falls back to filename inference
    // TODO: consume meta.seed once #615 wiring carries it (#588)
    // TODONT: database version injection. This can come out of versions.yml, but leave that to a cleaner database handling implementation.

    output:
    tuple val(meta), path("*.{mmcif,bcif}"), emit: modelcif
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // This will take a maximalist approach to populating the modelCIF
    // Every single structure from a sequence prediction method will be captured in .ModelGroup
    // Every common metric will be captured as a modelCIF qa_metric
    // The advantage is all software metadata and protocol steps are shared between all structures in the .ModelGroup
    // BinaryCIF output: pass --write_binary via task.ext.args.
    script:
    def args = task.ext.args ?: ''
    def container = task.container ?: 'None'
    """
    # The 'python:3.12-slim' base container only provides the interpreter; install the
    # pinned runtime dependencies declared in environment.yml so this works
    # identically under conda, docker and singularity. Skip the (re)install when
    # a conda environment has already provided them.
    python3 -c "import yaml, numpy, Bio.PDB, modelcif, msgpack" 2>/dev/null || \\
        pip install --quiet --no-cache-dir pyyaml==6.0.2 numpy biopython==1.84 "modelcif>=1.7" "msgpack>=1.0"

    populate_modelcif.py \\
        --structs ${structs} \\
        --msa ${msa} \\
        --plddt ${plddt} \\
        --pae ${pae} \\
        --ptm ${ptm} \\
        --iptm ${iptm} \\
        --name ${meta.id} \\
        --prog ${meta.model} \\
        --versions_yml ${versions_yml} \\
        --container_image ${container} \\
        --msa_tool ${meta.msa_tool ?: 'None'} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        modelcif: \$(python3 -c "import modelcif; print(modelcif.__version__)" 2>/dev/null || echo "unknown")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

}
