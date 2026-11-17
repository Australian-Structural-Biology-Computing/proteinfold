process MSA_CACHE_LOOKUP {
    tag "$meta.id"
    cpus 1
    memory '1 GB'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/183a5c5b5d21b977c9348a98e36d62610984a2017e427e6c182ef47782a0ea4f/data' :
        'community.wave.seqera.io/library/python:3.12--67983005212f7291' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path('cache_result')

    script:
    def cache_dir = MsaCache.directory(params) ?: ''
    def cache_id = params.msa_cache_id ?: ''
    def tool_id = params.msa_cache_tool_id ?: ''
    def args = task.ext.args ?: ''
    def artifact = task.ext.msa_cache_artifact ?: 'af3_json'
    def extension = artifact == 'a3m' ? 'a3m' : 'json'
    def quotedDir = "'${cache_dir.toString().replace("'", "'\\''")}'"
    def quotedId = "'${cache_id.toString().replace("'", "'\\''")}'"
    def quotedToolId = "'${tool_id.toString().replace("'", "'\\''")}'"
    def quotedArgs = "'${args.toString().replace("'", "'\\''")}'"
    def quotedFasta = "'${fasta.toString().replace("'", "'\\''")}'"
    """
    mkdir cache_result
    if [ -n ${quotedDir} ]; then
        msa_cache.py lookup ${quotedFasta} ${quotedId} ${quotedToolId} ${quotedArgs} ${quotedDir} cache_result/hit.${extension} ${artifact}
    fi
    """

    stub:
    def cache_dir = MsaCache.directory(params) ?: ''
    def cache_id = params.msa_cache_id ?: ''
    def tool_id = params.msa_cache_tool_id ?: ''
    def args = task.ext.args ?: ''
    def artifact = task.ext.msa_cache_artifact ?: 'af3_json'
    def extension = artifact == 'a3m' ? 'a3m' : 'json'
    def quotedDir = "'${cache_dir.toString().replace("'", "'\\''")}'"
    def quotedId = "'${cache_id.toString().replace("'", "'\\''")}'"
    def quotedToolId = "'${tool_id.toString().replace("'", "'\\''")}'"
    def quotedArgs = "'${args.toString().replace("'", "'\\''")}'"
    def quotedFasta = "'${fasta.toString().replace("'", "'\\''")}'"
    """
    mkdir cache_result
    if [ -n ${quotedDir} ]; then
        msa_cache.py lookup ${quotedFasta} ${quotedId} ${quotedToolId} ${quotedArgs} ${quotedDir} cache_result/hit.${extension} ${artifact}
    fi
    """
}
