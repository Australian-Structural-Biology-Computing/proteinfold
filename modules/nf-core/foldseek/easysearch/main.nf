process FOLDSEEK_EASYSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0':
        'biocontainers/foldseek:9.427df8a--pl5321hb365157_0' }"

    input:
    tuple val(meta)   , path(pdb)
    tuple val(meta_db), path(db)

    output:
    tuple val(meta), path("${meta.id}_${meta.mode_short ?: meta.model ?: meta.id}${meta.msa_short ? "_" + meta.msa_short : ""}.m8"), emit: aln, optional: true
    tuple val(meta), path("${meta.id}_${(meta.mode_short ?: meta.model ?: meta.id).toLowerCase()}${meta.msa_short ? "_" + meta.msa_short : ""}_foldseek.html"), emit: report, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def msa_tag = meta.msa_short ? "_${meta.msa_short}" : ''
    def output_file = "${meta.id}_${meta.mode_short ?: meta.model ?: meta.id}${msa_tag}.m8"
    if (args.contains("--format-mode 3")){
        output_file = "${meta.id}_${(meta.mode_short ?: meta.model ?: meta.id).toLowerCase()}${msa_tag}_foldseek.html"
    }
    """
    foldseek \\
        easy-search \\
        ${pdb} \\
        ${db}/${meta_db.id} \\
        ${output_file} \\
        tmpFolder \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.m8
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}
