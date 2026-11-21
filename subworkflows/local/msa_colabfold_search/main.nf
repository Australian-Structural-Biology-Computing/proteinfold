/*
 * MSA generation with colabfold_search (MMseqs2)
 *
 * Searches the samples against the ColabFold databases with MMseqs2, honouring
 * --msa_batch_size (samples per search task) and the MSA cache, and emits an
 * AlphaFold3-format MSA JSON per sample. AlphaFold3 JSON inputs are reduced to
 * their protein chains for the search and their original JSON is repopulated
 * with the generated MSAs afterwards (keeping PTMs, ligands, RNA and any other
 * non-protein entities intact).
 */

include { ALPHAFOLD3_JSON_TO_COLABFOLD_FASTA } from '../../../modules/local/alphafold3_json_to_colabfold_fasta'
include { MERGE_ALPHAFOLD3_MSA               } from '../../../modules/local/merge_alphafold3_msa'
include { MMSEQS_COLABFOLDSEARCH             } from '../../../modules/local/mmseqs_colabfoldsearch'
include { MMSEQS_COLABFOLDSEARCH_BATCH       } from '../../../modules/local/mmseqs_colabfoldsearch'
include { MULTIFASTA_TO_COLABFOLD_FASTA      } from '../../../modules/local/multifasta_to_colabfold_fasta'
include { MSA_CACHE_LOOKUP                   } from '../../../modules/local/msa_cache_lookup'

workflow MSA_COLABFOLD_SEARCH {

    take:
    ch_fasta         // channel: [meta, fasta] FASTA-origin samples
    ch_json          // channel: [meta, json] AlphaFold3-JSON-origin samples
    msa_batch_size   // int: number of samples processed together per search (-1 = one batch)
    ch_colabfold_db  // channel: path(colabfold_db)
    ch_uniref30      // channel: path(uniref30)

    main:
    // Mark JSON-origin samples so their search results are merged back into
    // the original (possibly rich) AlphaFold3 JSON afterwards.
    ch_json
        .map { meta, input_json ->
            def search_meta = meta.clone()
            search_meta.colabfold2_json_input = true
            [search_meta, input_json]
        }
        .set { ch_json_for_search }

    ch_json
        .map { meta, input_json -> [meta.id, input_json] }
        .set { ch_json_templates }

    ALPHAFOLD3_JSON_TO_COLABFOLD_FASTA(ch_json_for_search)

    // Colon-join multi-chain FASTAs so complexes are searched (and
    // paired) as a single query, matching the JSON extraction.
    MULTIFASTA_TO_COLABFOLD_FASTA(ch_fasta)

    ch_search_input = MULTIFASTA_TO_COLABFOLD_FASTA.out.query_fasta
        .mix(ALPHAFOLD3_JSON_TO_COLABFOLD_FASTA.out.query_fasta)

    def batch_size = msa_batch_size instanceof Number ? msa_batch_size : msa_batch_size.toString().toInteger()
    if (batch_size == 0 || batch_size < -1) {
        error("MSA batch size must be -1 or a positive integer, got ${batch_size}")
    }

    ch_generated_json = null
    def cache_dir = MsaCache.directory(params) ?: ''
    if (cache_dir) {
        // Cache lookups return the AF3-format MSA JSON directly.
        MSA_CACHE_LOOKUP(ch_search_input, cache_dir)
        MSA_CACHE_LOOKUP.out.branch { meta, csv, result ->
            hit: result.resolve('hit.json').exists()
            miss: true
        }.set { ch_cache_results }
        ch_search_input = ch_cache_results.miss.map { meta, csv, result -> [meta, csv] }
        ch_cached_json = ch_cache_results.hit.map { meta, csv, result -> [meta, result.resolve('hit.json')] }
    } else {
        ch_cached_json = channel.empty()
    }

    if (batch_size == 1) {
        MMSEQS_COLABFOLDSEARCH (
            ch_search_input,
            ch_colabfold_db,
            ch_uniref30,
            cache_dir
        )
        ch_generated_json = MMSEQS_COLABFOLDSEARCH.out.json
    } else {
        ch_batches = ch_search_input
            .map { meta, csv -> [meta.id, meta, csv] }
            .collect(flat: false)
            .flatMap { rows ->
                def sorted = rows.sort { a, b -> a[0] <=> b[0] }
                def ids = sorted.collect { row -> row[0] }
                if (ids.toSet().size() != ids.size()) error("Duplicate sample IDs are not allowed")
                def effective_batch_size = batch_size == -1 ? Math.max(1, sorted.size()) : batch_size
                sorted.collate(effective_batch_size).collect { batch ->
                    [batch.collect { item -> item[1] }, batch.collect { item -> item[2] }]
                }
            }
        MMSEQS_COLABFOLDSEARCH_BATCH(ch_batches, ch_colabfold_db, ch_uniref30, cache_dir)
        ch_generated_json = MMSEQS_COLABFOLDSEARCH_BATCH.out.json
            .flatMap { batch_meta, files ->
                def batch_files = files instanceof List ? files : [files]
                batch_files.collect { file ->
                    def sample_id = file.name.endsWith('.json') ? file.name[0..-6] : file.name
                    def sample = batch_meta.find { meta -> meta.get('id') == sample_id }
                    if (sample == null) error("MMseqs batch output ${file.name} does not match input samples ${batch_meta.collect { meta -> meta.id }}")
                    [sample, file]
                }
            }
    }

    // Cached hits and freshly generated results are treated uniformly.
    ch_search_result_all = ch_generated_json.mix(ch_cached_json)
    ch_search_result_all
        .branch { meta, msa_json ->
            json:  meta.colabfold2_json_input == true
            fasta: true
        }
        .set { ch_search_result }

    ch_search_result.json
        .map { meta, msa_json -> [meta.id, meta, msa_json] }
        .join(ch_json_templates)
        .map { id, meta, msa_json, original_json ->
            def clean_meta = meta.clone()
            clean_meta.remove('colabfold2_json_input')
            [clean_meta, msa_json, original_json]
        }
        .set { ch_json_merge_input }

    MERGE_ALPHAFOLD3_MSA(ch_json_merge_input)

    emit:
    msa_json = ch_search_result.fasta.mix(MERGE_ALPHAFOLD3_MSA.out.json)
}
