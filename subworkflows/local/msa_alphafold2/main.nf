/*
 * MSA generation with the alphafold2 data pipeline
 *
 * Runs alphafold's MSA-only pipeline per sample (the model preset is inferred
 * from the number of protein chains: single-chain entries run the monomer
 * pipeline, multi-chain entries the multimer pipeline with the uniprot
 * pairing search) and rebuilds an AlphaFold3-format MSA JSON with unpaired
 * and paired MSAs plus template hits. AlphaFold3 JSON inputs are reduced to
 * their protein chains and their original JSON is repopulated afterwards.
 */

include { ALPHAFOLD3_JSON_TO_AF2_FASTA } from '../../../modules/local/alphafold3_json_to_af2_fasta'
include { MERGE_ALPHAFOLD3_MSA         } from '../../../modules/local/merge_alphafold3_msa'
include { RUN_ALPHAFOLD2_MSA           } from '../../../modules/local/run_alphafold2_msa'
include { AF2_MSA_TO_COLABFOLD_JSON    } from '../../../modules/local/af2_msa_to_colabfold_json'
include { batchMergeInputs             } from '../utils_nfcore_proteinfold_pipeline'

def sanitise(id) {
    (id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
}


workflow MSA_ALPHAFOLD2 {

    take:
    ch_fasta                     // channel: [meta, fasta] FASTA-origin samples
    ch_json                      // channel: [meta, json] AlphaFold3-JSON-origin samples
    alphafold2_full_dbs          // boolean: full AF2 databases (otherwise reduced)
    alphafold2_max_template_date // string: latest template release date to keep
    uniref30_prefix              // string: prefix for uniref30 database files
    ch_af2_params                // channel: path(alphafold2 params)
    ch_af2_bfd                   // channel: path(bfd)
    ch_af2_small_bfd             // channel: path(small_bfd)
    ch_af2_mgnify                // channel: path(mgnify)
    ch_af2_pdb70                 // channel: path(pdb70)
    ch_af2_pdb_mmcif             // channel: path(pdb_mmcif)
    ch_af2_pdb_obsolete          // channel: path(pdb_obsolete)
    ch_af2_uniref30              // channel: path(uniref30)
    ch_af2_uniref90              // channel: path(uniref90)
    ch_af2_pdb_seqres            // channel: path(pdb_seqres)
    ch_af2_uniprot               // channel: path(uniprot)

    main:
    // AlphaFold3 JSON inputs: extract their protein chains into a
    // multi-record FASTA; the original JSON is repopulated with the generated
    // MSAs after the adaptor runs (keeping PTMs, ligands, RNA and any other
    // non-protein entities intact).
    ch_json
        .map { meta, input_json ->
            def json_meta = meta.clone()
            json_meta.af2_json_input = true
            [json_meta, input_json]
        }
        .set { ch_af2_json_input }

    ch_af2_json_by_id = ch_af2_json_input
        .map { meta, input_json -> [meta.id, meta, input_json] }

    ALPHAFOLD3_JSON_TO_AF2_FASTA(ch_af2_json_input)

    ch_af2_input = ch_fasta.mix(ALPHAFOLD3_JSON_TO_AF2_FASTA.out.query_fasta)

    RUN_ALPHAFOLD2_MSA (
        ch_af2_input,
        alphafold2_full_dbs,
        uniref30_prefix,
        ch_af2_params,
        ch_af2_bfd,
        ch_af2_small_bfd,
        ch_af2_mgnify,
        ch_af2_pdb70,
        ch_af2_pdb_mmcif,
        ch_af2_pdb_obsolete,
        ch_af2_uniref30,
        ch_af2_uniref90,
        ch_af2_pdb_seqres,
        ch_af2_uniprot
    )

    ch_af2_input
        .join(RUN_ALPHAFOLD2_MSA.out.raw_dir)
        .set { ch_af2_adaptor_input }
    AF2_MSA_TO_COLABFOLD_JSON(
        ch_af2_adaptor_input,
        !alphafold2_full_dbs,
        ch_af2_pdb_mmcif,
        ch_af2_pdb_obsolete,
        alphafold2_max_template_date
    )

    // FASTA-origin samples consume the adaptor JSON directly;
    // JSON-origin samples have their original JSON repopulated.
    AF2_MSA_TO_COLABFOLD_JSON.out.json
        .branch { meta, msa_json ->
            json:  meta.af2_json_input == true
            fasta: true
        }
        .set { ch_af2_adaptor_result }

    ch_af2_adaptor_result.json
        .map { meta, msa_json -> [meta.id, meta, msa_json] }
        .join(ch_af2_json_by_id)
        .map { id, meta, msa_json, json_meta, original_json ->
            def clean_meta = meta.clone()
            clean_meta.remove('af2_json_input')
            [clean_meta, msa_json, original_json]
        }
        .set { ch_af2_json_merge_input }

    // The alphafold2 path keeps its merge per sample (batch size 1) for now:
    // it runs through the same batched process, but its granularity stays
    // fixed and independent of --msa_batch_size and --analysis_batch_size.
    MERGE_ALPHAFOLD3_MSA(batchMergeInputs(ch_af2_json_merge_input, 1))
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
    msa_json = ch_af2_adaptor_result.fasta.mix(ch_merged_json)
}
