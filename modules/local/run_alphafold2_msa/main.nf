/*
 * Run Alphafold2 MSA
 *
 * Always runs the alphafold2 multimer data pipeline: the unpaired
 * uniref90/bfd/mgnify searches for every chain, the uniprot jackhmmer search
 * for heteromers (MSA pairing) and the pdb_seqres/hmmsearch template search.
 * Monomer and homomer entries skip the uniprot search inside the pipeline
 * itself.
 */
process RUN_ALPHAFOLD2_MSA {
    tag   "$meta.id"
    label 'process_medium'

    container "ghcr.io/jscgh/proteinfold_alphafold2_msa:2.1.0dev"

    input:
    tuple val(meta), path(fasta)
    val   db_preset
    val   uniref30_prefix
    path ('params/*')
    path ('bfd/*')
    path ('small_bfd/*')
    path ('mgnify/*')
    path ('pdb70/*')
    path ('pdb_mmcif/mmcif_files')
    path ('pdb_mmcif/*')
    path ('uniref30/*')
    path ('uniref90/*')
    path ('pdb_seqres/*')
    path ('uniprot/*')

    output:
    path ("raw/**")                             , emit: raw
    tuple val(meta), path ("raw")               , emit: raw_dir
    tuple val(meta), path ("raw/*/features.pkl"), emit: features
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>/dev/null | sed 's/Python //g' || echo \"unknown\""), emit: versions_python, topic: versions
    // The alphafold git repo is owned by the uid the image was built with, so
    // git's dubious-ownership check would fail when the container runs as any
    // other non-root uid (e.g. docker's '-u $(id -u)'). Exempt it via a global
    // config written under a throwaway HOME so no real HOME is touched.
    tuple val("${task.process}"), val('alphafold2'), eval("(export HOME=\$(mktemp -d); git config --global --add safe.directory /app/alphafold && cd /app/alphafold && git rev-parse HEAD) 2>/dev/null || echo \"unknown\""), emit: versions_alphafold2, topic: versions
    tuple val("${task.process}"), val('numpy'), eval("python3 -c \"import numpy; print(numpy.__version__)\" 2>/dev/null || echo \"unknown\""), emit: versions_numpy, topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python3 -c \"import Bio; print(Bio.__version__)\" 2>/dev/null || echo \"unknown\""), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ALPHAFOLD2_MSA module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def db_preset_cmd = db_preset ? "full_dbs --bfd_database_path=./bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --uniref30_database_path=./uniref30/${uniref30_prefix}" :
        "reduced_dbs --small_bfd_database_path=./small_bfd/bfd-first_non_consensus_sequences.fasta"
    """
    fix_obsolete.py pdb_mmcif/obsolete.dat > clean_obsolete.dat

    ## Handles multiple versions of mgnify database and selects the latest version
    mgnify_db_path=\$(ls -v ./mgnify/mgy_clusters*.fa | tail -n 1)

    ## The alphafold2 model preset is inferred per sample: entries with more
    ## than one protein chain run the multimer data pipeline (uniprot pairing
    ## search and pdb_seqres/hmmsearch templates), single-chain entries run the
    ## monomer pipeline (pdb70/hhsearch templates).
    n_chains=\$(grep -c '^>' "${fasta}")
    if [ "\$n_chains" -gt 1 ]; then
        model_preset=multimer
        extra_dbs="--pdb_seqres_database_path=./pdb_seqres/pdb_seqres.txt --uniprot_database_path=./uniprot/uniprot.fasta"
    else
        model_preset=monomer
        extra_dbs="--pdb70_database_path=./pdb70/pdb70"
    fi

    PYTHONPATH=/app/alphafold python3 run_msa.py \
        --fasta_paths=${fasta} \
        --model_preset=\$model_preset \$extra_dbs \
        --db_preset=${db_preset_cmd} \
        --output_dir=\$PWD \
        --data_dir=\$PWD \
        --uniref90_database_path=./uniref90/uniref90.fasta \
        --mgnify_database_path=\$mgnify_db_path \
        --template_mmcif_dir=./pdb_mmcif/mmcif_files \
        --obsolete_pdbs_path=./clean_obsolete.dat \
        $args

    # Can't use fasta.baseName to batch outputs in publishDir
    mv "${fasta.baseName}" raw/
    """

    stub:
    """
    mkdir -p ./raw/${fasta.baseName}/msas/A
    touch ./raw/${fasta.baseName}/features.pkl
    touch ./raw/${fasta.baseName}/msas/chain_id_map.json
    """
}
