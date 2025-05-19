process MMSEQS_CUDA {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    container "/home/z3545907/containers/mmseqs2.sif"

    input:
    tuple val(meta), path(fasta)
    path ('db/params')
    path colabfold_db
    path uniref30

    output:
    tuple val(meta), path("**.a3m"), emit: a3m
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local MMSEQS_CUDA module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def VERSION = '1.5.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ln -r -s $uniref30/uniref30_* ./db
    ln -r -s $colabfold_db/colabfold_envdb* ./db

     /usr/local/bin/mmseqs_avx2 easy-search \\
        ${fasta} \\
        ./db \\
        result.m8 \\
        tmp \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_search: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '1.5.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir results
    touch results/${meta.id}.a3m

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_search: $VERSION
    END_VERSIONS
    """
}
