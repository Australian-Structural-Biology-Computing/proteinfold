process COLABFOLD_BATCH {
    tag "$meta.id"
    label 'process_medium'

    container "nf-core/proteinfold_colabfold:dev"

    input:
    tuple val(meta), path(fasta)
    val   colabfold_model_preset
    path  ('params/*')
    path  ('colabfold_db/*')
    path  ('uniref30/*')
    val   numRec

    output:
    tuple val(meta), path ("${meta.id}_colabfold.pdb"), emit: top_ranked_pdb
    tuple val(meta), path ("*_relaxed_rank_*.pdb")    , emit: pdb
    tuple val(meta), path ("${meta.id}_plddt.tsv")    , emit: plddt
    tuple val(meta), path ("${meta.id}_msa.tsv")      , emit: msa
    tuple val(meta), path ("${meta.id}_*_pae.tsv")    , emit: paes
    tuple val(meta), path ("*_plddt.png")             , emit: plddt_img
    tuple val(meta), path ("*_coverage.png")          , emit: msa_img
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local COLABFOLD_BATCH module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def VERSION = '1.5.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ln -r -s params/alphafold_params_*/* params/
    colabfold_batch \\
        $args \\
        --num-recycle ${numRec} \\
        --data \$PWD \\
        --model-type ${colabfold_model_preset} \\
        ${fasta} \\
        \$PWD
    for i in `find *_relaxed_rank_001*.pdb`; do cp \$i `echo \$i | sed "s|_relaxed_rank_|\t|g" | cut -f1`"_colabfold.pdb"; done
    for i in `find *.png -maxdepth 0`; do cp \$i \${i%'.png'}_plddt.png; done
    cp *_relaxed_rank_001*.pdb ${meta.id}_colabfold.pdb

    extract_metrics.py  --name ${meta.id} \\
        --structs ${meta.id}_colabfold.pdb \\
        --a3ms  ${meta.id}/${meta.id}.a3m \\
        --jsons ${meta.id}/${meta.id}_scores_rank_*_alphafold2_ptm_model_*_seed_*.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_batch: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '1.5.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch "${meta.id}_colabfold.pdb"
    touch "${meta.id}_plddt.png"
    touch "${meta.id}_coverage.png"
    touch "${meta.id}_plddt.tsv"
    touch "${meta.id}_msa.tsv"
    touch "${meta.id}_0_pae.tsv"
    touch "${meta.id}_relaxed_rank_01.pdb"
    touch "${meta.id}_relaxed_rank_01.pdb"
    touch "${meta.id}_relaxed_rank_02.pdb"
    touch "${meta.id}_relaxed_rank_03.pdb"
    touch "${meta.id}_scores_rank.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_batch: $VERSION
    END_VERSIONS
    """
}
