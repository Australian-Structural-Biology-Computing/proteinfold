process GENERATE_REPORT {
    tag   "$meta.id-$meta.model"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/241f0746484727a3633f544c3747bfb77932e1c8c252e769640bd163232d9112/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pip_plotly:35975fa0fc54b2d3' }"

    input:
    tuple val(meta), path(pdb), path(msa, stageAs: 'input_msa'), path(pae, stageAs: 'input_pae'), path(iptm, stageAs: 'input_iptm'), path(ipsae, stageAs: 'input_ipsae'), path(chainwise, stageAs: 'input_chainwise'), path(chainwise_ipsae, stageAs: 'input_chainwise_ipsae')
    path(template)

    output:
    tuple val(meta), path ("*report.html")     , emit: report
    tuple val(meta), path ("*seq_coverage.png"), optional: true, emit: sequence_coverage
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('generate_report.py'), eval("python3 --version | sed 's/Python //g'"), emit: versions_generate_report, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    generate_report.py \\
        --type ${meta.mode_short} \
        --tag ${meta.mode_short}_${meta.msa_short ?: 'msa'} \\
        --msa ${msa} \\
        --pae ${pae} \\
        --iptm ${iptm} \\
        --ipsae ${ipsae} \\
        --chainwise ${chainwise} \\
        --chainwise_ipsae ${chainwise_ipsae} \\
        --pdb ${pdb.join(' ')} \\
        --html_template ${template} \\
        --output_dir ./ \\
        --name ${meta.id} \\
        $args \\
    """

    stub:
    """
    touch ${meta.id}_${meta.mode_short}_${meta.msa_short ?: 'msa'}_report.html
    touch ${meta.id}_${meta.mode_short}_${meta.msa_short ?: 'msa'}_seq_coverage.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_report.py: \$(python3 --version)
    END_VERSIONS
    """
}
