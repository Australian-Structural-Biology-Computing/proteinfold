/*
 * MSA generation with the alphafold3 data pipeline
 *
 * Runs alphafold3's data pipeline (MSA + template search, no inference) per
 * sample and emits the validated AlphaFold3-format MSA JSON.
 */

include { FASTA_TO_ALPHAFOLD3_JSON    } from '../../../modules/local/fasta_to_alphafold3_json'
include { RUN_ALPHAFOLD3_DATAPIPELINE } from '../../../modules/local/run_alphafold3_datapipeline'
include { VALIDATE_ALPHAFOLD3_MSA_JSON } from '../../../modules/local/validate_af3_msa_json'

workflow MSA_ALPHAFOLD3 {

    take:
    ch_fasta           // channel: [meta, fasta] FASTA-origin samples
    ch_json            // channel: [meta, json] AlphaFold3-JSON-origin samples
    ch_af3_small_bfd   // channel: path(small_bfd)
    ch_af3_mgnify      // channel: path(mgnify)
    ch_af3_mmcif_files // channel: path(mmcif_files)
    ch_af3_uniref90    // channel: path(uniref90)
    ch_af3_pdb_seqres  // channel: path(pdb_seqres)
    ch_af3_uniprot     // channel: path(uniprot)
    ch_af3_nt_rna      // channel: path(nt_rna)
    ch_af3_rfam        // channel: path(rfam)
    ch_af3_rnacentral  // channel: path(rnacentral)

    main:
    FASTA_TO_ALPHAFOLD3_JSON(ch_fasta)
    ch_af3_job_json = ch_json.mix(FASTA_TO_ALPHAFOLD3_JSON.out.json)

    RUN_ALPHAFOLD3_DATAPIPELINE (
        ch_af3_job_json,
        ch_af3_small_bfd,
        ch_af3_mgnify,
        ch_af3_mmcif_files,
        ch_af3_uniref90,
        ch_af3_pdb_seqres,
        ch_af3_uniprot,
        ch_af3_nt_rna,
        ch_af3_rfam,
        ch_af3_rnacentral
    )

    VALIDATE_ALPHAFOLD3_MSA_JSON(
        RUN_ALPHAFOLD3_DATAPIPELINE.out.data_json
    )

    emit:
    msa_json = VALIDATE_ALPHAFOLD3_MSA_JSON.out.json
}
