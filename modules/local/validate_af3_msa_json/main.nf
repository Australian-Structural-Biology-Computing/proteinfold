process VALIDATE_ALPHAFOLD3_MSA_JSON {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/biopython_python:06582b7b722f3db3"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_af3_msa_json.py \\
        ${json} \\
        --output ${meta.id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
