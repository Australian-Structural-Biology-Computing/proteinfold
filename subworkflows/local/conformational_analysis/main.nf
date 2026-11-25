//
// Conformational analysis of the predicted structures
//

include { FOLDSEEK_MULTIMERSEARCH } from '../../../modules/local/foldseek_multimersearch'
include { CLUSTER_CONFORMATIONS   } from '../../../modules/local/cluster_conformations'

workflow CONFORMATIONAL_ANALYSIS {

    take:
    ch_structures   // channel: tuple(meta, path structures) emitted per prediction mode
    ch_scores       // channel: tuple(meta, path scores JSONs) emitted per prediction mode
    cluster_cutoff  // float: complex TM-score similarity cutoff for conformational clustering

    main:
    // Every mode tags its outputs with its own meta map, so structures
    // predicted for the same sample arrive as separate per-mode tuples.
    // Flatten each tuple to one entry per structure file, tagged with the
    // mode that produced it, then regroup everything by sample id so the
    // analysis sees all predicted structures of a sample at once.
    // A glob output matching a single file arrives as a bare Path rather
    // than a list (and a Path iterates over its name components), so
    // normalise both to lists before flattening.
    ch_structures_by_sample = ch_structures
        .flatMap { meta, files -> (files instanceof List ? files : [files]).collect { file -> tuple(meta.id, [model: meta.model, mode_short: meta.mode_short], file) } }
        .groupTuple()

    ch_scores_by_sample = ch_scores
        .flatMap { meta, files -> (files instanceof List ? files : [files]).collect { file -> tuple(meta.id, file) } }
        .groupTuple()

    // The inference task errors out unless it produced both structures and
    // confidence scores, so the two grouped channels cover the same samples.
    ch_samples = ch_structures_by_sample
        .join(ch_scores_by_sample)
        .map { id, tags, structures, scores -> tuple([id: id], structures, scores, tags) }

    // Foldseek-Multimer aligns every pair of predicted structures of the
    // sample with an optimised chain pairing, giving chain-order-invariant
    // complex TM-scores (and the chain assignments used for each pair) for
    // clustering and 2D projection.
    FOLDSEEK_MULTIMERSEARCH(
        ch_samples.map { meta, structures, scores, tags -> tuple(meta, structures) }
    )

    CLUSTER_CONFORMATIONS(
        ch_samples
            .join(FOLDSEEK_MULTIMERSEARCH.out.multimer_report)
            .map { meta, structures, scores, tags, report -> tuple(meta, structures, scores, report, tags) },
        cluster_cutoff
    )

    emit:
    html            = CLUSTER_CONFORMATIONS.out.html
    table           = CLUSTER_CONFORMATIONS.out.table
    multimer_report = FOLDSEEK_MULTIMERSEARCH.out.multimer_report
}
