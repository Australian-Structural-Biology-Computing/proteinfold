process CIFCHECK {
    tag "$meta.id-$meta.model"
    label 'process_single'

    // NOTE: hard-depends on bin/dicts/mmcif_ma.sdb, a large binary ModelArchive
    // dictionary that is git-ignored and not bundled with this repository (see
    // .gitignore). Any invocation of this process requires that file to be
    // present on disk; it is expected to be supplied via nf-core/test-datasets
    // or a local setup step. The nf-test specs that depend on it are skipped
    // in nf-test.config until that data is wired in.
    input:
    tuple val(meta), path(mmcif)

    output:
    tuple val(meta), path(mmcif), emit: modelcif
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for mmcif_file in ${mmcif}; do
        CifCheck -f "\${mmcif_file}" -dictSdb ${projectDir}/bin/dicts/mmcif_ma.sdb
        if [ -s "\${mmcif_file}-diag.log" ] || [ -s "\${mmcif_file}-parser.log" ]; then
            echo "ModelArchive CifCheck validation errors in \${mmcif_file}:" >&2
            [ -s "\${mmcif_file}-diag.log" ]   && cat "\${mmcif_file}-diag.log"   >&2
            [ -s "\${mmcif_file}-parser.log" ] && cat "\${mmcif_file}-parser.log" >&2
            exit 1
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CifCheck: "2.500"
    END_VERSIONS
    """
}
