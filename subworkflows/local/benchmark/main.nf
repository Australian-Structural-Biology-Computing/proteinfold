//
// Reference-structure benchmarking for the top-ranked predicted structures
//

include { BENCHMARK_STRUCTURE } from '../../../modules/local/benchmark_structure'

workflow BENCHMARK {

    take:
    ch_predictions           // channel: [val(meta), path(structure)] of the top-ranked predicted models
    ch_references            // channel: [val(id), path(reference)] of the reference structures
    colabfold2_public_mode   // val: the selected ColabFold2 backend mode, or null
    benchmark_models         // val: normalized method names expected to emit predictions
    benchmark_strict         // val: abort the run when a referenced sample/method produces no benchmark

    main:
    // The pipeline config supplies a Boolean, but callers and tests may pass a CLI-style String, so
    // coerce an omitted or null value to false.
    def strict_benchmark = benchmark_strict instanceof Boolean
        ? benchmark_strict
        : Boolean.parseBoolean(benchmark_strict?.toString() ?: 'false')

    ch_normalized_predictions = ch_predictions.map { meta, structure ->
        def benchmark_meta = meta.clone()
        if (benchmark_meta.model == 'colabfold2') benchmark_meta.model = colabfold2_public_mode
        [benchmark_meta, structure]
    }

    // Keep only predictions for samples that declare a reference: combine pairs each reference
    // with every prediction of that sample and drops unreferenced samples, so they never reach
    // benchmarking and cannot trip its cardinality checks.
    ch_referenced_predictions = ch_normalized_predictions
        .map { meta, structure -> tuple(meta.id, meta, structure) }
        .combine(ch_references, by: 0)
        .map { _id, meta, structure, reference -> tuple(meta, structure, reference) }

    // Aggregate by (sample, method) so repeated emissions of the same pair collapse into a single
    // list and are detected here, rather than each launching its own BENCHMARK_STRUCTURE task and
    // silently colliding on the same published benchmark/<model>/<id>.json path.
    ch_grouped_predictions = ch_referenced_predictions
        .map { meta, structure, reference -> tuple([meta.id, meta.model], meta, [structure], [reference]) }
        .groupTuple()
        .map { key, metas, structures, references ->
            def first_meta = metas[0]
            if (metas.any { it != first_meta }) {
                error("Benchmarking received conflicting metadata for the same sample and method: ${key}")
            }
            def files = structures.flatten()
            if (files.size() != 1) {
                error("Benchmarking requires one predicted structure per sample and method: ${key[0]}, ${key[1]}")
            }
            tuple(first_meta, files[0], references[0])
        }

    BENCHMARK_STRUCTURE(ch_grouped_predictions)

    // Collate the per-sample JSON outputs into a deterministic, user-facing table so benchmark
    // scores can be compared across modes without parsing each OpenStructure result separately.
    ch_summary = BENCHMARK_STRUCTURE.out.scores
        .map { meta, score ->
            def result = new groovy.json.JsonSlurper().parse(score.toFile())
            tuple(meta.id, meta.model, result.lddt, result.tm_score)
        }
        .collect(flat: false)
        .filter { rows -> rows }
        .map { rows ->
            def csv_cell = { value ->
                def text = value?.toString() ?: ''
                text.contains(',') || text.contains('"') || text.contains('\n') || text.contains('\r')
                    ? '"' + text.replace('"', '""') + '"'
                    : text
            }
            def sorted_rows = rows.sort { left, right ->
                (left[0].toString() <=> right[0].toString()) ?: (left[1].toString() <=> right[1].toString())
            }
            def csv_rows = sorted_rows.collect { row -> row.collect(csv_cell).join(',') }
            ([ 'sample,mode,lddt,tm_score' ] + csv_rows).join('\n') + '\n'
        }
        .collectFile(
            storeDir: "${params.outdir}/benchmark",
            name: 'benchmark_summary.csv',
            newLine: false
        )

    ch_expected_benchmark_pairs = ch_references
        .combine(channel.fromList(benchmark_models.unique()))
        .map { id, _reference, model ->
            tuple([id, model], id, model)
        }
        .unique()

    ch_observed_benchmark_pairs = ch_grouped_predictions
        .map { meta, _structure, _reference ->
            tuple([meta.id, meta.model], true)
        }
        .unique()

    ch_missing_benchmark_pairs = ch_expected_benchmark_pairs
        // Unmatched expected pairs are released only after both channels close.
        .join(ch_observed_benchmark_pairs, by: 0, remainder: true)
        .filter { row ->
            row.size() == 4 && row[1] != null && row[3] == null
        }
        .map { _pair, id, model, _observed ->
            tuple(id, model)
        }

    if (!strict_benchmark) {
        ch_missing_benchmark_pairs.subscribe { id, model ->
            log.warn("Benchmark skipped: sample '${id}', method '${model}' produced no top-ranked structure.")
        }
    }

    // Persist the benchmark matrix gaps as durable output next to the published scores rather than
    // console warnings only, so an incomplete run can be detected from the results directory after
    // the fact. Nothing is written when every referenced prediction is present.
    ch_missing_report = ch_missing_benchmark_pairs
        .map { id, model -> "${id}\t${model}" }
        .collect()
        .filter { rows -> rows }
        .map { rows -> ([ 'sample\tmethod' ] + rows).join('\n') + '\n' }
        .collectFile(
            storeDir: "${params.outdir}/benchmark",
            name: 'benchmark_missing.tsv',
            newLine: false
        )

    if (strict_benchmark) {
        // collect() releases a single emission once all missing pairs are known, so the run aborts
        // deterministically with the complete gap list instead of racing per-pair errors.
        ch_missing_benchmark_pairs
            .map { id, model -> "${id}\t${model}" }
            .collect()
            .subscribe { rows ->
                if (rows) {
                    error(
                        'Strict benchmarking enabled: no top-ranked structure for the following referenced sample and method pairs: ' +
                        "${rows.join(', ')}. Set --benchmark_strict false to treat missing benchmarks as warnings."
                    )
                }
            }
    }

    emit:
    scores         = BENCHMARK_STRUCTURE.out.scores
    summary        = ch_summary
    missing        = ch_missing_benchmark_pairs
    missing_report = ch_missing_report
}
