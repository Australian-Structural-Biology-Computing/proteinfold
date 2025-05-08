/*
 * Run Alphafold2 PRED
 */
process RUN_ALPHAFOLD2_PRED {
    tag   "$meta.id"
    label 'process_medium'

    container "nf-core/proteinfold_alphafold2_split:dev"

    input:
    tuple val(meta), path(fasta)
    val   alphafold2_model_preset
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
    tuple val(meta2), path(features)

    output:
    path ("${fasta.baseName}*")
    tuple val(meta), path ("${meta.id}_alphafold2.pdb")              , emit: top_ranked_pdb
    tuple val(meta), path ("${fasta.baseName}/ranked*.pdb")          , emit: pdb
    tuple val(meta), path ("${fasta.baseName}/${meta.id}_plddt.tsv") , emit: plddt
    tuple val(meta), path ("${fasta.baseName}/${meta.id}_msa.tsv")   , emit: msa
    // TODO: alphafold2_model_preset == "monomer" the pae file won't exist, recommend running monomer_ptm by default. Performance loss tiny for insight
    // TODO: handle by passing 5 NO_FILE s
    tuple val(meta), path ("${fasta.baseName}/${meta.id}_*_pae.tsv") , emit: paes
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ALPHAFOLD2_PRED module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    """
    if [ -d params/alphafold_params_* ]; then ln -r -s params/alphafold_params_*/* params/; fi
    python3 /app/alphafold/run_predict.py \
        --fasta_paths=${fasta} \
        --model_preset=${alphafold2_model_preset} \
        --output_dir=\$PWD \
        --data_dir=\$PWD \
        --msa_path=${features} \
        $args

    cp "${fasta.baseName}"/ranked_0.pdb ./"${meta.id}"_alphafold2.pdb
    cd "${fasta.baseName}"

    extract_metrics.py --name ${meta.id} \\
        --structs *.pdb \\
        --pkls ${features} \\

    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ./"${meta.id}"_alphafold2.pdb
    touch ./"${meta.id}"_plddt.tsv
    touch ./${meta.id}_msa.tsv
    touch ./${meta.id}_0_pae.tsv
    mkdir "${fasta.baseName}"
    touch "${fasta.baseName}/ranked_0.pdb"
    touch "${fasta.baseName}/ranked_1.pdb"
    touch "${fasta.baseName}/ranked_2.pdb"
    touch "${fasta.baseName}/ranked_3.pdb"
    touch "${fasta.baseName}/ranked_4.pdb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
