//
// Download all the required chai1 parameters
//

include { ARIA2_UNCOMPRESS as ARIA2_CHAI1_CONFORMERS  } from './aria2_uncompress'
include { ARIA2_UNCOMPRESS as ARIA2_CHAI1_MODELS      } from './aria2_uncompress'
include { ARIA2_UNCOMPRESS as ARIA2_CHAI1_ESM         } from './aria2_uncompress'

workflow PREPARE_CHAI1_DBS {

    take:
    chai1_db                           // directory: /path/to/chai1/db/
    chai1_conformers_path              // directory: /path/to/db/
    chai1_models_path                  // directory: /path/to/db/models_v2/
    chai1_esm_path                     // directory: /path/to/db/esm/
    chai1_conformers_link              // directory: /link/to/db/
    chai1_models_link                  // directory: /link/to/db/models_v2/
    chai1_esm_link                     // directory: /link/to/db/esm/

    main:
    ch_versions   = Channel.empty()

    if (chai1_db) {
        ch_chai1_conformers   = Channel.value(files(chai1_conformers_path))
        ch_chai1_models       = Channel.value(files(chai1_models_path))
        ch_chai1_esm          = Channel.value(files(chai1_esm_path))
    }
    else {
        ARIA2_CHAI1_MODELS(
            chai1_models_link
        )
        ch_versions = ch_versions.mix(ARIA2_CHAI1_MODELS.out.versions)
        ch_chai1_models = ARIA2_CHAI1_MODELS.out.db

        ARIA2_CHAI1_CONFORMERS(
            chai1_conformers_link
        )
        ch_versions = ch_versions.mix(ARIA2_CHAI1_CONFORMERS.out.versions)
        ch_chai1_conformers = ARIA2_CHAI1_CONFORMERS.out.db

        ARIA2_CHAI1_ESM(
            chai1_esm_link
        )
        ch_versions = ch_versions.mix(ARIA2_CHAI1_ESM.out.versions)
        ch_chai1_esm = ARIA2_CHAI1_ESM.out.db

        ch_versions = ch_versions.mix(ARIA2_CHAI1_MODELS.out.versions)
    }

    emit:
    chai1_conformers    = ch_chai1_conformers
    chai1_models        = ch_chai1_models
    chai1_esm           = ch_chai1_esm
    versions            = ch_versions
}
