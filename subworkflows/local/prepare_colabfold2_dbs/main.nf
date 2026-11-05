include { MMSEQS_CREATEINDEX as MMSEQS_CREATEINDEX_COLABFOLDDB } from '../../../modules/nf-core/mmseqs/createindex/main'
include { MMSEQS_CREATEINDEX as MMSEQS_CREATEINDEX_UNIPROT30   } from '../../../modules/nf-core/mmseqs/createindex/main'
include { ARIA2_UNCOMPRESS as ARIA2_COLABFOLD_DB } from '../aria2_uncompress/main'
include { ARIA2_UNCOMPRESS as ARIA2_UNIREF30     } from '../aria2_uncompress/main'

process DOWNLOAD_COLABFOLD2_WEIGHTS {
    tag "$colabfold2_model ($colabfold2_precision)"
    label 'process_medium'
    container "ghcr.io/jscgh/proteinfold_colabfold:af3-preview-e9e8bd0"

    input:
    val colabfold2_model
    val colabfold2_precision

    output:
    path 'colabfold2_models', emit: params_path

    script:
    def markerRevision = colabfold2_model == 'chai1' ? '_r3' : (colabfold2_model == 'boltz2' ? '_r1' : '')
    """
    set -euo pipefail
    mkdir -p colabfold2_models
    python3 - <<'PY'
    from pathlib import Path
    from colabfold.alphafold3.models import ensure_weights
    from colabfold.alphafold3.weights import ensure_ccd
    root = Path('colabfold2_models')
    ensure_ccd(root)
    ensure_weights('${colabfold2_model}', data_dir=root, precision='${colabfold2_precision}', download=True)
    PY
    test -f colabfold2_models/params/af3/libcifpp/components.cif
    test -f colabfold2_models/params/af3/${colabfold2_model}${colabfold2_precision == 'int8' ? '-int8' : ''}/download_${colabfold2_model}_${colabfold2_precision}${markerRevision}_finished.txt
    """
}

workflow PREPARE_COLABFOLD2_DBS {
    take:
    colabfold2_params_path
    colabfold_db
    colabfold2_model
    colabfold2_precision
    use_msa_server
    colabfold_envdb_path
    colabfold_uniref30_path
    colabfold_db_link
    colabfold_uniref30_link
    colabfold_create_index

    main:
    def localModels = colabfold_db && !colabfold_db.toString().contains('://') ? "${colabfold_db}/colabfold2_models" : null
    def companions = [
        'esmfold2_lm300m': 'esmc_300m',
        'esmfold2_lm600m': 'esmc_600m'
    ]
    def suffix = colabfold2_precision == 'int8' ? '-int8' : ''
    def companion = companions[colabfold2_model]
    def markerRevision = colabfold2_model == 'chai1' ? '_r3' : (colabfold2_model == 'boltz2' ? '_r1' : '')
    def marker = "params/af3/${colabfold2_model}${suffix}/download_${colabfold2_model}_${colabfold2_precision}${markerRevision}_finished.txt"
    def localReady = localModels && new File(localModels, marker).isFile() && new File(localModels, 'params/af3/libcifpp/components.cif').isFile()
    def rosettafold3Blob = 'params/af3/rosettafold3/rosettafold3.bin.zst'
    def rosettafold3Bytes = 1361361731L
    if (localReady && colabfold2_model == 'rosettafold3' && colabfold2_precision == 'fp32') {
        localReady = new File(localModels, rosettafold3Blob).length() == rosettafold3Bytes
    }
    if (localReady && companion) {
        localReady = new File(localModels, "params/af3/${colabfold2_model}${suffix}/${colabfold2_model}.lm.npz").isFile() && new File(localModels, "params/af3/${companion}/${companion}.bin.zst").isFile()
    }
    def weightsRoot = colabfold2_params_path ?: (localReady ? localModels : null)
    def af3_dir = "${weightsRoot}/params/af3"
    def model_dir = "${af3_dir}/${colabfold2_model}${suffix}"

    if (weightsRoot) {
        file(weightsRoot, type: 'dir', checkIfExists: true)
        if (!workflow.stubRun) {
            file(model_dir, type: 'dir', checkIfExists: true)
            file("${model_dir}/download_${colabfold2_model}_${colabfold2_precision}${markerRevision}_finished.txt", checkIfExists: true)
            if (colabfold2_model == 'rosettafold3' && colabfold2_precision == 'fp32' && !weightsRoot.toString().contains('://') && new File(weightsRoot.toString(), rosettafold3Blob).length() != rosettafold3Bytes) {
                error "RosettaFold3 weights in ${weightsRoot} do not match the checkpoint required by this ColabFold2 image (${rosettafold3Bytes} bytes). Use a refreshed weights cache."
            }
            if (companion) {
                file("${model_dir}/${colabfold2_model}.lm.npz", checkIfExists: true)
                file("${af3_dir}/${companion}/${companion}.bin.zst", checkIfExists: true)
            }
            file("${af3_dir}/libcifpp/components.cif", checkIfExists: true)
        }
        ch_params_path = channel.value(file(weightsRoot, type: 'dir', checkIfExists: true))
    } else if (workflow.stubRun) {
        ch_params_path = channel.value(file('assets/dummy_db_dir', type: 'dir', checkIfExists: true))
    } else {
        DOWNLOAD_COLABFOLD2_WEIGHTS(colabfold2_model, colabfold2_precision)
        ch_params_path = DOWNLOAD_COLABFOLD2_WEIGHTS.out.params_path
    }

    ch_colabfold_db = channel.empty()
    ch_uniref30 = channel.empty()
    if (!use_msa_server) {
        if (colabfold_db) {
            ch_colabfold_db = channel.value(files(colabfold_envdb_path, type: 'any', checkIfExists: true))
            ch_uniref30 = channel.value(files(colabfold_uniref30_path, type: 'any', checkIfExists: true))
        } else {
            ARIA2_COLABFOLD_DB(colabfold_db_link)
            ch_colabfold_db = ARIA2_COLABFOLD_DB.out.db

            if (colabfold_create_index) {
                MMSEQS_CREATEINDEX_COLABFOLDDB (
                    ch_colabfold_db.map { path_str ->
                        def db_file = file(path_str)
                        [[id: 'colabfolddb'], db_file]
                    }
                )
                ch_colabfold_db = MMSEQS_CREATEINDEX_COLABFOLDDB.out.db_indexed.map { _meta, dir -> files("${dir}/*") }
            } else {
                ch_colabfold_db = ch_colabfold_db.map { dir_path -> files("${dir_path}/*") }
            }

            ARIA2_UNIREF30(colabfold_uniref30_link)
            ch_uniref30 = ARIA2_UNIREF30.out.db

            if (colabfold_create_index) {
                MMSEQS_CREATEINDEX_UNIPROT30 (
                    ch_uniref30.map { path_str ->
                        def db_file = file(path_str)
                        [[id: 'uniprot30'], db_file]
                    }
                )
                ch_uniref30 = MMSEQS_CREATEINDEX_UNIPROT30.out.db_indexed.map { _meta, dir -> files("${dir}/*") }
            } else {
                ch_uniref30 = ch_uniref30.map { dir_path -> files("${dir_path}/*") }
            }
        }
    }

    emit:
    params_path = ch_params_path
    colabfold_db = ch_colabfold_db
    uniref30 = ch_uniref30
}
