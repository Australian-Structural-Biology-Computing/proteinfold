process GENERATE_MULTIQC_CONTENTS {
    tag   "$meta.model"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/183a5c5b5d21b977c9348a98e36d62610984a2017e427e6c182ef47782a0ea4f/data' :
        'community.wave.seqera.io/library/python:3.12--67983005212f7291' }"

    input:
    tuple val(meta), path(metric_files)
    path(generator_script)

    output:
    tuple val(meta), path("*_mqc.json"), emit: mqc_json
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('generate_multiqc_contents.py'), eval("python3 --version | sed 's/Python //g'"), emit: versions_generator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${generator_script} \\
        --model ${meta.model} \\
        --output-dir ./ \\
        ${metric_files.join(' ')} \\
        $args
    """

    stub:
    """
    touch proteinfold_${meta.model}_generalstats_mqc.json proteinfold_${meta.model}_plddt_lineplot_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_multiqc_contents.py: \$(python3 --version)
    END_VERSIONS
    """
}
