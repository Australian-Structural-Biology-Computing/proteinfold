process TO_COLABFOLD_FASTA {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    tuple val(meta), path(input_file, stageAs: 'input/*')
    val cache_dir // absolute MSA cache directory; '' disables the lookup

    output:
    tuple val(meta), path("${meta.id}.fasta"), path("cache_result"), emit: query

    when:
    task.ext.when == null || task.ext.when

    script:
    // The input is staged under its own name inside input/, keeping it
    // apart from the <id>.fasta output even when a FASTA input is itself
    // named after the sample. The cache lookup keys on the query content
    // and must hash the same arguments as the store side of the search
    // (--db-load-mode), otherwise cached searches are never hit.
    def staged = input_file instanceof Collection ? input_file[0].name : input_file.name
    def input_name = staged.tokenize('/')[-1]
    def dir = (cache_dir ?: '').toString().replace("'", "'\\''")
    def cache_id = (params.msa_cache_id ?: '').toString().replace("'", "'\\''")
    def tool_id = (params.msa_cache_tool_id ?: '').toString().replace("'", "'\\''")
    def cache_args = params.colabfold_db_load_mode ? "--db-load-mode ${params.colabfold_db_load_mode}" : ''
    def cache_args_quoted = cache_args.toString().replace("'", "'\\''")
    """
    to_colabfold_query_fasta.py \\
        input/${input_name} \\
        --id ${meta.id} \\
        --output ${meta.id}.fasta

    mkdir cache_result
    if [ -n '${dir}' ]; then
        msa_cache.py lookup '${meta.id}.fasta' '${cache_id}' '${tool_id}' '${cache_args_quoted}' '${dir}' cache_result/hit.json af3_json
    fi
    """

    stub:
    // Stub runs must produce the real colon-joined query and run the real
    // lookup: downstream cache hits key on the query content even in
    // stub-mode test runs, where a placeholder would turn every hit into a
    // miss. Both helpers only need the host's standard library.
    def staged = input_file instanceof Collection ? input_file[0].name : input_file.name
    def input_name = staged.tokenize('/')[-1]
    def dir = (cache_dir ?: '').toString().replace("'", "'\\''")
    def cache_id = (params.msa_cache_id ?: '').toString().replace("'", "'\\''")
    def tool_id = (params.msa_cache_tool_id ?: '').toString().replace("'", "'\\''")
    def cache_args = params.colabfold_db_load_mode ? "--db-load-mode ${params.colabfold_db_load_mode}" : ''
    def cache_args_quoted = cache_args.toString().replace("'", "'\\''")
    """
    to_colabfold_query_fasta.py \\
        input/${input_name} \\
        --id ${meta.id} \\
        --output ${meta.id}.fasta

    mkdir cache_result
    if [ -n '${dir}' ]; then
        msa_cache.py lookup '${meta.id}.fasta' '${cache_id}' '${tool_id}' '${cache_args_quoted}' '${dir}' cache_result/hit.json af3_json
    fi
    """
}
