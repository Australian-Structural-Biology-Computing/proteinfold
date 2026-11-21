process ALPHAFOLD3_JSON_TO_COLABFOLD_A3M {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    tuple val(meta), path(input_json)

    output:
    tuple val(meta), path("${meta.id}.a3m"), path("*compatibility.flag"), emit: a3m
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if af3_msa_json_to_a3m.py \\
        ${input_json} \\
        --output ${meta.id}.a3m; then
        touch compatible.compatibility.flag
    else
        conversion_status=\$?
        if [ "\$conversion_status" -ne 2 ] && [ "\$conversion_status" -ne 3 ]; then
            echo "af3_msa_json_to_a3m.py failed with status \$conversion_status" >&2
            exit "\$conversion_status"
        fi
        # Exit 2 means the JSON has no reusable protein MSA (no protein
        # entity, or no unpairedMsa): emit an empty a3m and mark the sample
        # incompatible with the AlphaFold2 backends instead of failing the
        # run, so multi-mode invocations can still fold it elsewhere.
        touch ${meta.id}.a3m
        touch incompatible.compatibility.flag
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    if af3_msa_json_to_a3m.py \\
        ${input_json} \\
        --output ${meta.id}.a3m; then
        touch compatible.compatibility.flag
    else
        conversion_status=\$?
        if [ "\$conversion_status" -ne 2 ] && [ "\$conversion_status" -ne 3 ]; then
            exit "\$conversion_status"
        fi
        touch ${meta.id}.a3m
        touch incompatible.compatibility.flag
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
