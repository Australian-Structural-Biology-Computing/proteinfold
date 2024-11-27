//
// TBD: Download all the required Rosettafold-All-Atom databases and parameters
//


workflow PREPARE_ROSETTAFOLD_ALL_ATOM_DBS {

    take:
    bfd_rosettafold_all_atom_path      // directory: /path/to/bfd/
    uniref30_rosettafold_all_atom_path // directory: /path/to/uniref30/rosettafold_all_atom/
    pdb100_rosettafold_all_atom_path
    rfaa_paper_weights_path

    main:
    ch_bfd                  = Channel.value(file(bfd_rosettafold_all_atom_path))
    ch_uniref30             = Channel.value(file(uniref30_rosettafold_all_atom_path))
    ch_pdb100               = Channel.value(file(pdb100_rosettafold_all_atom_path))
    ch_rfaa_paper_weights   = Channel.value(file(rfaa_paper_weights_path))
    ch_versions             = Channel.empty()

    emit:
    bfd                 = ch_bfd
    uniref30            = ch_uniref30
    pdb100              = ch_pdb100
    rfaa_paper_weights  = ch_rfaa_paper_weights
    versions            = ch_versions
}
