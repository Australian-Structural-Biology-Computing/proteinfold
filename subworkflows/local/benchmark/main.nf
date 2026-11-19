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

    main:
    ch_normalized_predictions = ch_predictions.map { meta, structure ->
        def benchmark_meta = meta.clone()
        if (benchmark_meta.model == 'colabfold2') benchmark_meta.model = colabfold2_public_mode
        [benchmark_meta, structure]
    }

    ch_benchmark_predictions = ch_normalized_predictions.map { meta, structure ->
        [meta.id, meta, structure]
    }
    BENCHMARK_STRUCTURE(
        ch_benchmark_predictions
            .combine(ch_references, by: 0)
            .map { id, meta, structures, reference ->
                def files = structures instanceof List ? structures : [structures]
                if (files.size() != 1) {
                    error("Benchmarking requires one predicted structure per sample and method: ${id}, ${meta.model}")
                }
                [meta, files[0], reference]
            }
    )

    ch_expected_benchmark_pairs = ch_references
        .combine(channel.fromList(benchmark_models.unique()))
        .map { id, _reference, model ->
            tuple([id, model], id, model)
        }
        .unique()

    ch_observed_benchmark_pairs = ch_normalized_predictions
        .map { meta, _structure ->
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

    ch_missing_benchmark_pairs.subscribe { id, model ->
        log.warn("Benchmark skipped: sample '${id}', method '${model}' produced no top-ranked structure.")
    }

    emit:
    scores = BENCHMARK_STRUCTURE.out.scores
    missing = ch_missing_benchmark_pairs
}
