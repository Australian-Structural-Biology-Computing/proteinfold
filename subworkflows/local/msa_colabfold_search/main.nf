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

    // One conversion task converts every sample into its ColabFold query
    // FASTA and consults the MSA cache in the same task. It emits the miss
    // queries pre-batched: each chunk FASTA concatenates msa_batch_size
    // sample records and feeds exactly one search task, and every cache hit
    // leaves only its cached AF3-format MSA JSON behind.
    def batch_size = msa_batch_size instanceof Number ? msa_batch_size : msa_batch_size.toString().toInteger()
    if (batch_size == 0 || batch_size < -1) {
        error("MSA batch size must be -1 or a positive integer, got ${batch_size}")
    }
    def cache_dir = MsaCache.directory(params) ?: ''

    TO_COLABFOLD_FASTA(
        ch_fasta.mix(ch_json_for_search)
            .collect(flat: false)
            .filter { rows -> !rows.isEmpty() }
            // Several samples may reference the same input file; stage it
            // once and let the task map every sample onto its staged name.
            .map { rows -> [rows, rows.collect { it[1] }.unique { it.toString() }] },
        batch_size,
        cache_dir,
        false
    )

    // Cache hits: match each cached AF3-format MSA JSON back to its sample
    // by the sanitised-id prefix of the file name.
    ch_cached_json = TO_COLABFOLD_FASTA.out.hits
        .flatMap { pairs, files ->
            // Private copies: emitted lists may share buffers across
            // concurrent consumers.
            def list = new ArrayList(files instanceof List ? files : [files])
            list.collect { file ->
                def sid = file.name.replaceAll(/\.json$/, '')
                def pair = pairs.find { (it[0].id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_') == sid }
                pair == null ? error("Cached MSA result ${file.name} does not match any sample in the batch") : [pair[0], file]
            }
        }

    // One search task per chunk FASTA; the sample metadata list rides along
    // so each result can be matched back to its sample by record header.
    MMSEQS_COLABFOLDSEARCH_BATCH(
        TO_COLABFOLD_FASTA.out.misses
            .flatMap { pairs, files ->
                def list = new ArrayList(files instanceof List ? files : [files])
                list.collect { file -> [new ArrayList(pairs.collect { it[0] }), file] }
            },
        ch_colabfold_db,
        ch_uniref30,
        cache_dir
    )
    ch_generated_json = MMSEQS_COLABFOLDSEARCH_BATCH.out.json
        .flatMap { samples, files ->
            def list = new ArrayList(files instanceof List ? files : [files])
            list.collect { file ->
                def sid = file.name.endsWith('.json') ? file.name[0..-6] : file.name
                def sample = samples.find { (it.id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_') == sid }
                if (sample == null) error("MMseqs search result ${file.name} does not match any input sample")
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
