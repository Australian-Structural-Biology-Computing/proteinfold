process CLUSTER_CONFORMATIONS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/241f0746484727a3633f544c3747bfb77932e1c8c252e769640bd163232d9112/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pip_plotly:35975fa0fc54b2d3' }"

    input:
    tuple val(meta), path(structures), path(scores), path(multimer_report), val(model_labels)
    val cluster_cutoff

    output:
    tuple val(meta), path ("*conformations.html"), emit: html
    tuple val(meta), path ("*conformations.tsv") , emit: table
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('cluster_conformations.py'), eval("python3 --version | sed 's/Python //g'"), emit: versions_cluster_conformations, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // One manifest row per structure so the script knows which prediction
    // model (public mode) produced it; the lists are index-aligned by
    // construction in the calling subworkflow.
    def manifest = model_labels.withIndex().collect { label, i ->
        "${structures[i].name}\t${label.model}\t${label.mode_short}"
    }.join('\n')
    def scores_arg = scores instanceof List ? scores.join(' ') : scores
    def argv = [
        'cluster_conformations.py',
        "--name ${meta.id}",
        "--cluster-cutoff ${cluster_cutoff}",
        '--manifest manifest.tsv',
        "--multimer-report ${multimer_report}",
        "--structures ${structures instanceof List ? structures.join(' ') : structures}",
        scores_arg ? "--scores ${scores_arg}" : null,
        '--output-dir .',
        args,
    ].findAll { it }.join(' ')

    // The manifest is interpolated as a multi-line value, which disables
    // Nextflow's script-block dedent; the heredoc terminators therefore stay
    // at column 0 so bash recognises them regardless of dedent behaviour.
    """
    cat > manifest.tsv <<MANIFEST
${manifest}
MANIFEST

    ${argv}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //g')
    cluster_conformations.py: \$(python3 --version | sed 's/Python //g')
END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_conformations.html
    touch ${meta.id}_conformations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        cluster_conformations.py: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
