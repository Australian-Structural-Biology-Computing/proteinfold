/*
 * Run Alphafold2 MSA
 */
process RUN_ALPHAFOLD2_MSA {
    tag   "$meta.id"
    label 'process_medium'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ALPHAFOLD2_MSA module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    container "nf-core/proteinfold_alphafold2_msa:dev"

    input:
    tuple val(meta), path(fasta)
    val   db_preset
    val   alphafold2_model_preset
    path ('params/*')
    path ('bfd/*')
    path ('small_bfd/*')
    path ('mgnify/*')
    path ('pdb70/*')
    path ('pdb_mmcif/*')
    path ('uniref30/*')
    path ('uniref90/*')
    path ('pdb_seqres/*')
    path ('uniprot/*')

    output:
    path ("${fasta.baseName}*")
    tuple val(meta), path ("${fasta.baseName}.features.pkl"), emit: features
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_preset = db_preset ? "full_dbs --bfd_database_path=${params.alphafold2_db}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --uniref30_database_path=${params.alphafold2_db}/uniref30/UniRef30_2021_03" :
        "reduced_dbs --small_bfd_database_path=${params.alphafold2_db}/small_bfd/bfd-first_non_consensus_sequences.fasta"
    if (alphafold2_model_preset == 'multimer') {
        alphafold2_model_preset += " --pdb_seqres_database_path=${params.alphafold2_db}/pdb_seqres/pdb_seqres.txt --uniprot_database_path=${params.alphafold2_db}/uniprot/uniprot.fasta "
    }
    else {
        alphafold2_model_preset += " --pdb70_database_path=${params.alphafold2_db}/pdb70/pdb70_from_mmcif_200916/pdb70 "
    }
    """
    RUNTIME_TMP=\$(mktemp -d)
    if [ -f ${params.alphafold2_db}/pdb_seqres/pdb_seqres.txt ]
        cp ${params.alphafold2_db}/pdb_seqres/pdb_seqres.txt \${RUNTIME_TMP}
        then sed -i "/^\\w*0/d" \${RUNTIME_TMP}/pdb_seqres.txt
    fi
    python3 /app/alphafold/run_msa.py \
        --fasta_paths=${fasta} \
        --model_preset=${alphafold2_model_preset} \
        --db_preset=${db_preset} \
        --output_dir=\$PWD \
        --data_dir=\$PWD \
        --uniref90_database_path=${params.alphafold2_db}/uniref90/uniref90.fasta \
        --mgnify_database_path=${params.alphafold2_db}/mgnify/mgy_clusters_2022_05.fa \
        --template_mmcif_dir=${params.alphafold2_db}/pdb_mmcif/mmcif_files \
        --obsolete_pdbs_path=${params.alphafold2_db}/pdb_mmcif/obsolete.dat  \
        $args

    cp "${fasta.baseName}"/features.pkl ./"${fasta.baseName}".features.pkl
    rm -rf "\${RUNTIME_TMP}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ./"${fasta.baseName}".features.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(gawk --version| head -1 | sed 's/GNU Awk //; s/, API:.*//')
    END_VERSIONS
    """
}
