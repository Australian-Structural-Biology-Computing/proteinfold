process FOLDSEEK_MULTIMERSEARCH {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/88bcddb476149bd710cbc6b61b11ca5ab724efe78c0b6977813de92e9ffc160d/data' :
        'community.wave.seqera.io/library/foldseek:10.941cd33--519bd079ce8d3fc7' }"

    input:
    tuple val(meta), path(structures, stageAs: 'structures/*')

    output:
    tuple val(meta), path("${meta.id}_multimer_report.tsv"), emit: multimer_report
    tuple val("${task.process}"), val('foldseek'), eval("foldseek --help |& sed -n 's/.*Version: //p'"), emit: versions_foldseek, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // All-vs-all complex alignment over every predicted structure of the
    // sample. Foldseek-Multimer optimises the chain pairing between two
    // complexes, so the resulting TM-scores are invariant to chain naming,
    // ordering and composition differences between prediction backends, and
    // the report records the assignment it aligned each pair with.
    // Exhaustive search skips the prefilter: sample-level result sets are
    // small (models x seeds x backends) and a complete distance matrix is
    // required for clustering and projection.
    """
    foldseek \\
        easy-multimersearch \\
        structures/ \\
        structures/ \\
        multimer_result \\
        tmp \\
        --exhaustive-search 1 \\
        --chain-name-mode 1 \\
        --threads ${task.cpus} \\
        $args

    mv multimer_result_report "${meta.id}_multimer_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help |& sed -n 's/.*Version: //p')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}_multimer_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: 941cd33ff0771cd2e3f144e3293e22a2b87e9fda
    END_VERSIONS
    """
}
