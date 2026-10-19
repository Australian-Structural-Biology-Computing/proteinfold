process CIFCHECK {
    tag "$meta.id-$meta.model"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/keiran-rowell-unsw/proteinfold-modelcif:1.7'

    input:
    tuple val(meta), path(mmcif), path(ma_dictionary)

    output:
    tuple val(meta), path(mmcif), emit: modelcif
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for modelcif_file in ${mmcif}; do
        case "\${modelcif_file}" in
            *.bcif)
                ${projectDir}/bin/modelcif_bcif_to_mmcif.py "\${modelcif_file}" "\${modelcif_file}.cifcheck.mmcif"
                check_file="\${modelcif_file}.cifcheck.mmcif"
                ;;
            *) check_file="\${modelcif_file}" ;;
        esac
        ${projectDir}/bin/CifCheck -f "\${check_file}" -dictSdb ${ma_dictionary}
        if [ -s "\${check_file}-diag.log" ] || [ -s "\${check_file}-parser.log" ]; then
            cat "\${check_file}-diag.log" "\${check_file}-parser.log" >&2
            exit 1
        fi
    done
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        CifCheck: "2.500"
    END_VERSIONS
    """
}
