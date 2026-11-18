//
// Reference-structure benchmarking for the top-ranked predicted structures
//

include { BENCHMARK_STRUCTURE } from '../../../modules/local/benchmark_structure'

workflow BENCHMARK {

    take:
    ch_predictions           // channel: [val(meta), path(structure)] of the top-ranked predicted models
    ch_references            // channel: [val(id), path(reference)] of the reference structures
    colabfold2_public_mode   // val: the selected ColabFold2 backend mode, or null

    main:
    ch_benchmark_predictions = ch_predictions.map { meta, structure ->
        [meta.id, meta, structure]
    }
    BENCHMARK_STRUCTURE(
        ch_benchmark_predictions
            .combine(ch_references, by: 0)
            // KNOWN LIMITATION: combine silently drops references that have no
            // matching prediction (e.g. an optional top-ranked output that never
            // emits, as in ColabFold2). Such samples get no benchmark and the
            // pipeline still succeeds. TODO: reconcile expected (id, model) pairs
            // against emitted scores and fail on any shortfall.
            .map { id, meta, structures, reference ->
                def files = structures instanceof List ? structures : [structures]
                if (files.size() != 1) {
                    error("Benchmarking requires one predicted structure per sample and method: ${id}, ${meta.model}")
                }
                def benchmark_meta = meta.clone()
                if (benchmark_meta.model == 'colabfold2') benchmark_meta.model = colabfold2_public_mode
                [benchmark_meta, files[0], reference]
            }
    )

    emit:
    scores = BENCHMARK_STRUCTURE.out.scores
}
