process AF2_MSA_TO_COLABFOLD_JSON {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/jscgh/proteinfold_alphafold2_msa:2.1.0"

    input:
    tuple val(meta), path(fasta), path(raw_dir)
    val use_small_bfd

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def small_bfd_arg = use_small_bfd ? '--small-bfd' : ''
    """
    af2_msa_to_colabfold_json.py \\
        --fasta ${fasta} \\
        --msa_dir raw_dir \\
        --output ${meta.id}.json \\
        --name ${meta.id} \\
        ${small_bfd_arg} \\
        $args

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
