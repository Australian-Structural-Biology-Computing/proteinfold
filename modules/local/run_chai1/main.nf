process RUN_CHAI1 {
    tag "$meta.id"
    label 'process_medium'

    container "jscrh/chai:dev"

    input:
    tuple val(meta), path(fasta)
    path('*')
    path('models_v2/*')
    path('esm')

    output:
    tuple val(meta), path ("${meta.id}_chai1.cif")    , emit: pdb
    tuple val(meta), path ("${meta.id}_plddt_mqc.tsv"), emit: multiqc
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_CHAI1 module does not support Conda. Please use Docker / Singularity / Podman / Apptainer instead.")
    }
    def args = task.ext.args ?: ''
    """
    export CHAI_DOWNLOADS_DIR="./"
    mamba run -n chai1 chai-lab fold \
        ${meta.id}.fasta \
        ${meta.id} \
        $args

    cp ${meta.id}/pred.model_idx_0.cif ./${meta.id}_chai1.cif
    cd ${meta.id}
    awk '{print \$6"\\t"\$11}' pred.model_idx_0.cif > ranked_0_plddt.tsv
    for i in 1 2 3 4
        do awk '{print \$6"\\t"\$11}' pred.model_idx_\$i.cif | awk '{print \$2}' > ranked_"\$i"_plddt.tsv
    done
    paste ranked_0_plddt.tsv ranked_1_plddt.tsv ranked_2_plddt.tsv ranked_3_plddt.tsv ranked_4_plddt.tsv > plddt.tsv
    echo -e Positions"\\t"rank_0"\\t"rank_1"\\t"rank_2"\\t"rank_3"\\t"rank_4 > header.tsv
    cat header.tsv plddt.tsv > ../"${meta.id}"_plddt_mqc.tsv
    cd ../

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        chai_lab: \$(conda run -n chai1 python -c "import chai_lab; print(chai_lab.__version__)")
        torch: \$(conda run -n chai1 python -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ./"${meta.id}"_chai1.cif
    touch ./"${meta.id}"_plddt_mqc.tsv
    mkdir ./"${meta.id}"
    touch "${meta.id}/ranked_0.pdb"
    touch "${meta.id}/ranked_1.pdb"
    touch "${meta.id}/ranked_2.pdb"
    touch "${meta.id}/ranked_3.pdb"
    touch "${meta.id}/ranked_4.pdb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        chai_lab: \$(conda run -n chai1 python -c "import chai_lab; print(chai_lab.__version__)")
        torch: \$(conda run -n chai1 python -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """
}
