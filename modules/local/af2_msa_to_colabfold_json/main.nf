/*
 * Convert the alphafold2 MSA + template pipeline output into an
 * AlphaFold3-format input JSON
 *
 * Rebuilds the colabfold_batch AF3 job JSON from RUN_ALPHAFOLD2_MSA: per-chain
 * unpaired MSAs (uniref90 -> bfd -> mgnify), the pipeline_multimer species
 * pairing for heteromers (uniref90 reuse for homomers), and - when the
 * template mmCIFs are available - the chain template hits as alphafold3
 * folding_input.Template entries. Runs in the colabfold container because the
 * template conversion reuses colabfold's build_templates (alphafold3).
 */
process AF2_MSA_TO_COLABFOLD_JSON {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    input:
    tuple val(meta), path(fasta), path(raw_dir, stageAs: 'raw_dir')
    val   use_small_bfd
    path(template_mmcif_dir, stageAs: 'template_mmcif_dir')
    path(obsolete_dat, stageAs: 'obsolete.dat')
    val max_template_date

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local AF2_MSA_TO_COLABFOLD_JSON module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def small_bfd_arg = use_small_bfd ? '--small-bfd' : ''
    """
    af2_msa_to_colabfold_json.py \\
        --fasta ${fasta} \\
        --msa_dir raw_dir \\
        --output ${meta.id}.json \\
        --name ${meta.id} \\
        --template_mmcif_dir template_mmcif_dir \\
        --obsolete obsolete.dat \\
        --max-template-date ${max_template_date} \\
        ${small_bfd_arg} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
