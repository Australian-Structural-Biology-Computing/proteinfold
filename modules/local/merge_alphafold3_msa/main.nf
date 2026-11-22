include { batchTag } from '../../../subworkflows/local/utils_nfcore_proteinfold_pipeline'

process MERGE_ALPHAFOLD3_MSA {
    tag   "${batchTag(batch_meta)}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    tuple val(batch_meta), path(msa_jsons, stageAs: 'searched/*'), path(original_jsons, stageAs: 'original/*'), val(manifest)
    // manifest: one '<sample id>\t<searched MSA file>\t<original input file>' line per sample

    output:
    tuple val(batch_meta), path("*.json"), emit: json
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // One task merges a whole batch: each row repopulates the sample's
    // original (possibly rich) AlphaFold3 JSON with its searched MSAs, and
    // the merged JSON is named after the sanitised sample id.
    def safe = { name -> "'" + (name as String).replaceAll("'", "'\\''") + "'" }
    """
cat > manifest.tsv <<'MERGE_MANIFEST'
${manifest}
MERGE_MANIFEST

while IFS=\$'\\t' read -r sid msa_name original_name <&3; do
    merge_alphafold3_msa.py \\
        "searched/\${msa_name}" \\
        "original/\${original_name}" \\
        --output "\${sid}.json"
done 3< manifest.tsv

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
"""

    stub:
    // The merge is pure-python JSON repopulation; the original input is a
    // real staged file even in stub runs, so hand it through unchanged to
    // keep downstream consumers (a3m conversion, inference) honest.
    """
cat > manifest.tsv <<'MERGE_MANIFEST'
${manifest}
MERGE_MANIFEST

while IFS=\$'\\t' read -r sid msa_name original_name <&3; do
    cp "original/\${original_name}" "\${sid}.json"
done 3< manifest.tsv

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
"""
}
