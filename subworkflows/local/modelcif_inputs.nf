//
// Build the ASSEMBLE_MODELCIF input tuple per sample+model (#615).
//
// Shape matches modules/local/assemble_modelcif:
//   [ meta, structs, msa, plddt, pae, ptm, iptm, versions ]
//
// Structures arrive as [meta, path] or [meta, [paths]] and are grouped by
// (id, model) into one rank-ordered list; each optional metric stream is
// joined on (id, model) and filled from the caller's dummy channel when a
// program does not produce it (e.g. ESMFold has no MSA/PAE/PTM/ipTM).
//
workflow MODELcif_INPUTS {

    take:
    ch_structs       // channel: [ meta, path(structs) ] (top-ranked or all ranks)
    ch_msa           // channel: [ meta, path(msa.tsv) ]
    ch_plddt         // channel: [ meta, path(plddt_mqc.tsv) ]
    ch_pae           // channel: [ meta, path(*_pae.tsv) ]
    ch_ptm           // channel: [ meta, path(*_ptm.tsv) ]
    ch_iptm          // channel: [ meta, path(*_iptm.tsv) ]
    ch_dummy_file    // channel: value(file(NO_FILE)) or [file] single value
    ch_versions_rows // channel: [ [ id:'*', model:'*', versions:{...} ] ]

    main:
    ch_keyed = ch_structs
        .map { m, p -> [ m.id, m.model, (p instanceof List) ? p : [ p ] ] }
        .groupTuple(by: [0, 1])

    ch_m   = ch_msa.map   { a, b, c -> [ a, b, c ] }
    ch_p1 = ch_plddt.map { a, b, c -> [ a, b, c ] }
    ch_p2 = ch_pae.map   { a, b, c -> [ a, b, c ] }
    ch_p3 = ch_ptm.map   { a, b, c -> [ a, b, c ] }
    ch_p4 = ch_iptm.map  { a, b, c -> [ a, b, c ] }

    ch_with_msa  = ch_keyed.join(ch_m,   by: [0, 1], failOnMismatch: false)
    ch_with_p1   = ch_with_msa.join(ch_p1, by: [0, 1], failOnMismatch: false)
    ch_with_p2   = ch_with_p1.join(ch_p2, by: [0, 1], failOnMismatch: false)
    ch_with_p3   = ch_with_p2.join(ch_p3, by: [0, 1], failOnMismatch: false)
    ch_with_p4   = ch_with_p3.join(ch_p4, by: [0, 1], failOnMismatch: false)

    ch_all = ch_with_p4
        .join(
            ch_versions_rows.map { [ '*', '*', it.versions ] },
            by: [0, 1],
            failOnMismatch: false,
        )

    ch_modelcif = ch_all.map { row ->
        def vals = row.drop(2).collect { x -> (x instanceof List && x.size() == 1) ? x[0] : x }
        [ [ id: row[0], model: row[1] ] ] + vals
    }

    emit:
    modelcif = ch_modelcif  // channel: [ meta, structs, msa, plddt, pae, ptm, iptm, versions ]
}
