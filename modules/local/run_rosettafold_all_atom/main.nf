/*
 * Run RoseTTAFold_All_Atom
 */
process RUN_ROSETTAFOLD_ALL_ATOM {
    tag "$meta.id"
    label 'process_medium'

    container "nf-core/proteinfold_rosettafold_all_atom:dev"

    input:
    tuple val(meta), path(fasta)
    path ('bfd/*')
    path ('UniRef30_2020_06/*')
    path ('pdb100_2021Mar03/*')
    path ('*')

    output:
    tuple val(meta), path ("${meta.id}_rosettafold_all_atom.pdb"), emit: pdb
    tuple val(meta), path ("${meta.id}_plddt.tsv")               , emit: plddt
    tuple val(meta), path ("${meta.id}_msa.tsv")                 , emit: msa
    tuple val(meta), path ("${meta.id}_*_pae.tsv")               , emit: paes
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ROSETTAFOLD_ALL_ATOM module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def VERSION = '1.2.0dev' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mamba run --name RFAA python /app/RoseTTAFold-All-Atom/rf2aa/run_inference.py \
    --config-dir /app/RoseTTAFold-All-Atom/rf2aa/config/inference \
    --config-name "${fasta}" \
    $args


    cp "${fasta.baseName}".pdb ./"${meta.id}"_rosettafold_all_atom.pdb

    extract_metrics.py --name ${meta.id} \\
        --structs "${meta.id}_rosettafold_all_atom.pdb" \\
        --a3ms "${fasta.baseName}/A/t000_.msa0.a3m" \\
        --pts ${meta.id}_aux.pt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}_rosettafold_all_atom.pdb"
    touch "${meta.id}.pdb"
    touch "${meta.id}_aux.pt"
    touch "${meta.id}_plddt.tsv"
    touch "${meta.id}_msa.tsv"
    touch "${meta.id}_0_pae.tsv"
    mkdir "${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
