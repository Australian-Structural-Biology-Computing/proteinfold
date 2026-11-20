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
    def quoted_model     = "'${model.toString().replace("'", "'\\''")}'"
    def quoted_reference = "'${reference.toString().replace("'", "'\\''")}'"
    """
    benchmark_structure.py \\
        --model ${quoted_model} \\
        --reference ${quoted_reference} \\
        --output benchmark.json
    """

    stub:
    """
    echo '{"lddt": 1.0, "tm_score": 1.0}' > benchmark.json
    """
}
