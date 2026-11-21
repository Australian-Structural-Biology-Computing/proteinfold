process TO_COLABFOLD_FASTA {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    tuple val(meta), path(input_file, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: query_fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    // The input is staged under its own name inside input/, keeping it
    // apart from the <id>.fasta output even when a FASTA input is itself
    // named after the sample.
    def staged = input_file instanceof Collection ? input_file[0].name : input_file.name
    def input_name = staged.tokenize('/')[-1]
    """
    to_colabfold_query_fasta.py \\
        input/${input_name} \\
        --id ${meta.id} \\
        --output ${meta.id}.fasta
    """

    stub:
    // Stub runs must produce the real colon-joined query: downstream cache
    // lookups key on the query content even in stub-mode test runs, where an
    // empty placeholder would turn every cache hit into a miss. The stub
    // therefore runs the same converter on the host (it only needs the
    // standard library).
    def staged = input_file instanceof Collection ? input_file[0].name : input_file.name
    def input_name = staged.tokenize('/')[-1]
    """
    to_colabfold_query_fasta.py \\
        input/${input_name} \\
        --id ${meta.id} \\
        --output ${meta.id}.fasta
    """
}
