process TO_COLABFOLD_FASTA {
    tag   "queries_${pairs.size()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    tuple val(pairs), path(input_files, stageAs: 'input/*')
    val msa_batch_size  // records per emitted chunk FASTA; -1 = one chunk with everything
    val emit_queries    // keep one query FASTA per sample (for consumers that need per-sample files)

    output:
    // Every chunk is a plain concatenation of msa_batch_size ColabFold
    // query FASTAs (one record per sample), ready for the search.
    tuple val(pairs), path("chunks/*.fasta"), optional: true, emit: chunks
    // One ColabFold query FASTA per sample, named after the sanitised sample
    // id; the record header inside carries the raw id. Only emitted when
    // emit_queries is true (e.g. the MSA server path folds one query per
    // sample).
    tuple val(pairs), path("queries/*.fasta"), optional: true, emit: queries

    when:
    task.ext.when == null || task.ext.when

    script:
    // pairs: one [meta, input file] entry per sample; input_files: the
    // deduplicated files behind them, so several samples may share one
    // staged input without a name collision.
    def metas = pairs.collect { it[0] }
    def files = input_files instanceof Collection ? input_files.collect() : [input_files]
    def sanitise = { id -> (id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    // sanitised sample id, raw record header, staged input name; one
    // tab-separated triple per line
    def manifest = pairs.collect { meta, file ->
        "${sanitise(meta.id)}\t${meta.id}\t${file.name}"
    }
    if (metas.collect(sanitise).toSet().size() != metas.size()) {
        error("TO_COLABFOLD_FASTA received sample IDs that collide after sanitisation: ${metas*.id}")
    }
    if (files.collect { it.name }.toSet().size() != files.size()) {
        error("TO_COLABFOLD_FASTA received distinct input files that share a name; rename them so one task can stage them together")
    }
    if ([metas*.id, files*.name].flatten().any { (it as String).matches('(?s).*[\\t\\r\\n].*') }) {
        error("TO_COLABFOLD_FASTA received sample IDs or input file names containing tabs or newlines")
    }
    def limit = msa_batch_size instanceof Number ? msa_batch_size : msa_batch_size.toString().toInteger()
    def mode_arg = emit_queries ? '--queries-dir queries' : '--chunks-dir chunks'
    def command = [
        'to_colabfold_query_fasta.py batch \\',
        '    --manifest manifest.tsv \\',
        "    --batch-size ${limit} \\",
        "    ${mode_arg}",
    ]

    // One batch-native call: converts every sample in memory and
    // concatenates the queries into pre-batched chunk FASTAs. No
    // per-sample intermediate files are written.
    """
cat > manifest.tsv <<'MANIFEST'
${manifest.join('\n')}
MANIFEST

${command.join('\n')}
rmdir chunks queries 2>/dev/null || true
"""
}

def sanitise(id) {
    (id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
}
