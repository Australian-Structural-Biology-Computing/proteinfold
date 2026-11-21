include { MMSEQS_CREATEINDEX as MMSEQS_CREATEINDEX_COLABFOLDDB } from '../../../modules/nf-core/mmseqs/createindex/main'
include { MMSEQS_CREATEINDEX as MMSEQS_CREATEINDEX_UNIPROT30   } from '../../../modules/nf-core/mmseqs/createindex/main'
include { ARIA2_UNCOMPRESS as ARIA2_COLABFOLD_DB } from '../aria2_uncompress/main'
include { ARIA2_UNCOMPRESS as ARIA2_UNIREF30     } from '../aria2_uncompress/main'

process DOWNLOAD_COLABFOLD2_WEIGHTS {
    tag "$colabfold2_model ($colabfold2_precision)"
    label 'process_medium'
    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    input:
    tuple val(colabfold2_model), val(colabfold2_precision)

    output:
    tuple val(colabfold2_model), path('colabfold2_models'), emit: params_path

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

process DOWNLOAD_COLABFOLD2_AF2_WEIGHTS {
    tag "alphafold2 (${colabfold2_variants.join(',')})"
    label 'process_medium'
    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    input:
    val colabfold2_variants

    output:
    path 'af2_params', emit: params_path

    script:
    def variant_list = colabfold2_variants.collect { variant -> "'${variant}'" }.join(', ')
    """
    set -euo pipefail
    python3 - <<'PY'
    from pathlib import Path
    from colabfold.download import download_alphafold_params
    for variant in [${variant_list}]:
        download_alphafold_params(variant, data_dir=Path('af2_params'))
    PY
    test -d af2_params/params
    """

    stub:
    """
    mkdir -p af2_params/params
    touch af2_params/params/download_finished.txt
    """
}

// ColabFold's AlphaFold2 parameter cache: <root>/params holds the npz
// files plus the success markers its own downloader writes and skips on.
def af2CacheReady(root, variant) {
    def markers = [
        'alphafold2_ptm'        : ['download_finished.txt', 'alphafold_params'],
        'alphafold2'            : ['download_finished.txt', 'alphafold_params'],
        'alphafold2_multimer_v3': ['download_complexes_multimer_v3_finished.txt', 'alphafold_multimer_params'],
        'alphafold2_multimer_v2': ['download_complexes_multimer_v2_finished.txt', 'alphafold_multimer_params'],
        'alphafold2_multimer_v1': ['download_complexes_multimer_v1_finished.txt', 'alphafold_multimer_params'],
    ]
    def entry = markers[variant]
    if (entry == null) {
        return false
    }
    def params_dir = new File(root.toString(), 'params')
    def marker = new File(params_dir, entry[0])
    def npz_present = params_dir.listFiles()?.any { it.name.startsWith(entry[1]) && it.name.endsWith('.npz') }
    return marker.isFile() && npz_present
}

workflow PREPARE_COLABFOLD2_DBS {
    take:
    colabfold2_params_path
    colabfold_db
    colabfold2_models
    colabfold2_precision
    use_msa_server
    colabfold_envdb_path
    colabfold_uniref30_path
    colabfold_db_link
    colabfold_uniref30_link
    colabfold_create_index

    main:
    // Resolve the weights root once per requested model type and emit
    // [model_type, params_dir] tuples so the dispatch loop can pair each
    // mode with the weights it needs (AF2-family and AF3-family layouts
    // differ, and multiple AF3-family modes each need their own weights).
    def localModels = colabfold_db && !colabfold_db.toString().contains('://') ? "${colabfold_db}/colabfold2_models" : null
    def companions = [
        'esmfold2_lm300m': 'esmc_300m',
        'esmfold2_lm600m': 'esmc_600m'
    ]
    def suffix = colabfold2_precision == 'int8' ? '-int8' : ''
    def rosettafold3Blob = 'params/af3/rosettafold3/rosettafold3.bin.zst'
    def rosettafold3Bytes = 1361361731L

    def af2Variants = null
    def af2Models = []
    def af3Models = []
    colabfold2_models.each { model ->
        if (model == 'auto') {
            af2Variants = (af2Variants ?: []) + ['alphafold2_ptm', 'alphafold2_multimer_v3']
            af2Models << model
        } else if (model.startsWith('alphafold2')) {
            af2Variants = (af2Variants ?: []) + [model]
            af2Models << model
        } else {
            af3Models << model
        }
    }

    def ch_params_specs = channel.empty()

    // AF2-family weights use ColabFold's own layout (<root>/params with
    // alphafold_params .npz files); one root serves every AF2 variant, and
    // 'auto' needs both the ptm and the multimer_v3 parameter sets.
    if (af2Models) {
        if (colabfold2_params_path) {
            file(colabfold2_params_path, type: 'dir', checkIfExists: true)
            if (!workflow.stubRun) {
                file("${colabfold2_params_path}/params", type: 'dir', checkIfExists: true)
            }
            af2Models.each { model ->
                ch_params_specs = ch_params_specs.mix(channel.value([model, file(colabfold2_params_path, type: 'dir', checkIfExists: true)]))
            }
        } else if (workflow.stubRun) {
            af2Models.each { model ->
                ch_params_specs = ch_params_specs.mix(channel.value([model, file('assets/dummy_db_dir', type: 'dir', checkIfExists: true)]))
            }
        } else if (localModels && af2Variants.every { variant -> af2CacheReady(localModels, variant) }) {
            // Reuse a complete implicit cache: colabfold's own success
            // markers (the same ones its downloader skips on) plus the
            // parameter files gate the reuse, exactly like the AF3-family
            // branch above.
            file("${localModels}/params", type: 'dir', checkIfExists: true)
            af2Models.each { model ->
                ch_params_specs = ch_params_specs.mix(channel.value([model, file(localModels, type: 'dir', checkIfExists: true)]))
            }
        } else {
            DOWNLOAD_COLABFOLD2_AF2_WEIGHTS(af2Variants.unique())
            ch_params_specs = ch_params_specs.mix(
                DOWNLOAD_COLABFOLD2_AF2_WEIGHTS.out.params_path.map { dir -> af2Models.collect { [model: it, dir: dir] } }.flatten().map { spec -> [spec.model, spec.dir] }
            )
        }
    }

    // AF3-family weights live under <root>/params/af3/<model>; every
    // requested AF3-family model must be present in (or downloaded to) the
    // shared root.
    if (af3Models) {
        def weightsRoot = colabfold2_params_path ?: (localModels ? localModels : null)
        if (colabfold2_params_path || (localModels && af3Models.every { model ->
            def markerRevision = model == 'chai1' ? '_r3' : (model == 'boltz2' ? '_r1' : '')
            def marker = "params/af3/${model}${suffix}/download_${model}_${colabfold2_precision}${markerRevision}_finished.txt"
            def ready = new File(localModels, marker).isFile() && new File(localModels, 'params/af3/libcifpp/components.cif').isFile()
            if (ready && model == 'rosettafold3' && colabfold2_precision == 'fp32') {
                ready = new File(localModels, rosettafold3Blob).length() == rosettafold3Bytes
            }
            def companion = companions[model]
            if (ready && companion) {
                ready = new File(localModels, "params/af3/${model}${suffix}/${model}.lm.npz").isFile() && new File(localModels, "params/af3/${companion}/${companion}.bin.zst").isFile()
            }
            ready
        })) {
            def root = weightsRoot
            if (!workflow.stubRun) {
                af3Models.each { model ->
                    def markerRevision = model == 'chai1' ? '_r3' : (model == 'boltz2' ? '_r1' : '')
                    def model_dir = "${root}/params/af3/${model}${suffix}"
                    file(model_dir, type: 'dir', checkIfExists: true)
                    file("${model_dir}/download_${model}_${colabfold2_precision}${markerRevision}_finished.txt", checkIfExists: true)
                    if (model == 'rosettafold3' && colabfold2_precision == 'fp32' && !root.toString().contains('://') && new File(root.toString(), rosettafold3Blob).length() != rosettafold3Bytes) {
                        error "RosettaFold3 weights in ${root} do not match the checkpoint required by this ColabFold2 image (${rosettafold3Bytes} bytes). Use a refreshed weights cache."
                    }
                    def companion = companions[model]
                    if (companion) {
                        file("${model_dir}/${model}.lm.npz", checkIfExists: true)
                        file("${root}/params/af3/${companion}/${companion}.bin.zst", checkIfExists: true)
                    }
                }
                file("${root}/params/af3/libcifpp/components.cif", checkIfExists: true)
            }
            af3Models.each { model ->
                ch_params_specs = ch_params_specs.mix(channel.value([model, file(root, type: 'dir', checkIfExists: true)]))
            }
        } else if (workflow.stubRun) {
            af3Models.each { model ->
                ch_params_specs = ch_params_specs.mix(channel.value([model, file('assets/dummy_db_dir', type: 'dir', checkIfExists: true)]))
            }
        } else {
            // Channel-fed invocation: processes cannot be called from within
            // a closure, so each AF3-family model streams in as a tuple.
            ch_af3_weights = channel.fromList(af3Models).map { model -> tuple(model, colabfold2_precision) }
            DOWNLOAD_COLABFOLD2_WEIGHTS(ch_af3_weights)
            ch_params_specs = ch_params_specs.mix(DOWNLOAD_COLABFOLD2_WEIGHTS.out.params_path)
        }
    }

    ch_params_path = ch_params_specs

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
