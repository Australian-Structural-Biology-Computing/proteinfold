process MODELCIF_VALIDATE {
    tag "$meta.id-$meta.model"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "python:3.12-slim"

    input:
    tuple val(meta), path(mmcif)

    output:
    tuple val(meta), path(mmcif), emit: modelcif
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # The 'python:3.12-slim' base container only provides the interpreter; install the
    # pinned runtime dependency declared in environment.yml so this works
    # identically under conda, docker and singularity. Skip the (re)install when
    # a conda environment has already provided it.
    python3 -c "import modelcif" 2>/dev/null || \\
        pip install --quiet --no-cache-dir "modelcif>=1.7"

    python3 <<'PYEOF'
import warnings
import modelcif.reader
import sys

files = "${mmcif}".split()
for f in files:
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        with open(f) as fh:
            systems = modelcif.reader.read(fh)
    if not systems:
        raise ValueError(f"No ModelCIF data blocks found in {f}")
    for system in systems:
        if not system.entities:
            raise ValueError(f"ModelCIF system in {f} has no entities")
        if not system.protocols:
            raise ValueError(f"ModelCIF system in {f} has no modeling protocol (missing _ma_protocol_step)")
        if not system.model_groups:
            raise ValueError(f"ModelCIF system in {f} has no model groups (missing _ma_model_group / _ma_model_list)")
    print(f'py-modelcif validation passed: {f}', file=sys.stderr)

with open('versions.yml', 'w') as fh:
    import modelcif
    fh.write('${task.process}:\\n')
    fh.write(f'    modelcif: {modelcif.__version__}\\n')
    import platform
    fh.write(f'    python: {platform.python_version()}\\n')
PYEOF
    """

}
