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

include { TO_COLABFOLD_FASTA } from '../../../modules/local/to_colabfold_fasta'
include { MERGE_ALPHAFOLD3_MSA               } from '../../../modules/local/merge_alphafold3_msa'
include { batchMergeInputs                   } from '../utils_nfcore_proteinfold_pipeline'
include { MMSEQS_COLABFOLDSEARCH_BATCH       } from '../../../modules/local/mmseqs_colabfoldsearch'
include { MSA_CACHE_LOOKUP                   } from '../../../modules/local/msa_cache_lookup'

def sanitise(id) {
    (id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
}


workflow MSA_COLABFOLD_SEARCH {

    take:
    ch_fasta               // channel: [meta, fasta] FASTA-origin samples
    ch_json                // channel: [meta, json] AlphaFold3-JSON-origin samples
    msa_batch_size         // int: number of samples processed together per search (-1 = one batch)
    ch_colabfold_db        // channel: path(colabfold_db)
    ch_uniref30            // channel: path(uniref30)

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

    // Colon-join every sample into a single ColabFold query so complexes
    // are searched (and paired) as one entry: multi-chain FASTAs join
    // across their records, and AlphaFold3 JSONs reduce to their protein
    // chains, matching the JSON extraction.
    TO_COLABFOLD_FASTA(ch_fasta.mix(ch_json_for_search))

    ch_search_input = TO_COLABFOLD_FASTA.out.query_fasta

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

    // One batched search path serves every batch size: 1 runs one sample per
    // task, -1 runs all samples of the run in a single task. Each input is a
    // colon-delimited FASTA holding exactly one record whose header is the
    // sample id, so a batch is a plain concatenation and colabfold_search
    // names every result after its record header.
    ch_batches = ch_search_input
        .map { meta, csv -> [meta.id, meta, csv] }
        .collect(flat: false)
        .flatMap { rows ->
            // Private copies: the collected rows travel through nextflow
            // operators that may share buffers across concurrent consumers.
            def sorted = new ArrayList(rows).sort { a, b -> a[0] <=> b[0] }
            def ids = sorted.collect { row -> row[0] }
            if (ids.toSet().size() != ids.size()) error("Duplicate sample IDs are not allowed")
            def effective_batch_size = batch_size == -1 ? Math.max(1, sorted.size()) : batch_size
            sorted.collate(effective_batch_size).collect { batch ->
                [new ArrayList(batch.collect { item -> item[1] }), new ArrayList(batch.collect { item -> item[2] })]
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

    // One merge task per MSA batch: every row repopulates its sample's
    // original AlphaFold3 JSON with the searched MSAs, and the merged JSONs
    // are matched back to their samples by the sanitised-id file name. The
    // merge belongs to the MSA stage, so it groups by --msa_batch_size like
    // the search itself.
    MERGE_ALPHAFOLD3_MSA(batchMergeInputs(ch_json_merge_input, batch_size))
    ch_merged_json = MERGE_ALPHAFOLD3_MSA.out.json
        .flatMap { batch_meta, files ->
            def list = new ArrayList(files instanceof List ? files : [files])
            list.collect { file ->
                def sid = file.name.replaceAll(/\.json$/, '')
                def meta = batch_meta.find { sanitise(it.id) == sid }
                meta == null ? error("Merged MSA JSON ${file.name} does not match any sample in the batch") : [meta, file]
            }
        }

    emit:
    msa_json = ch_search_result.fasta.mix(ch_merged_json)
}
