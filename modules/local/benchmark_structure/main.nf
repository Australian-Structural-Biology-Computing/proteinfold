process BENCHMARK_STRUCTURE {
    tag "$meta.id-$meta.model"
    label 'process_single'

    container 'registry.scicore.unibas.ch/schwede/openstructure@sha256:9f2669ffde16c45939fc6e8e6dd4f907697301c081edf3d97d464c034a421d19'

    input:
    tuple val(meta), path(model, stageAs: 'model/*'), path(reference, stageAs: 'reference/*')

    output:
    tuple val(meta), path('benchmark.json'), emit: scores
    tuple val("${task.process}"), val('openstructure'), eval("ost --version | sed 's/OpenStructure //'"), emit: versions, topic: versions

    script:
    def quote = { value -> "'${value.toString().replace("'", "'\\''")}'" }
    """
    benchmark_structure.py \\
        --model ${quote(model)} \\
        --reference ${quote(reference)} \\
        --output benchmark.json
    """

    stub:
    """
    echo '{}' > benchmark.json
    """
}
