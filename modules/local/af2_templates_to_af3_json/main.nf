/*
 * Add alphafold2 template hits to an AlphaFold3-format MSA JSON
 *
 * Turns the pdb_seqres/hmmsearch template hits produced by RUN_ALPHAFOLD2_MSA
 * into alphafold3 folding_input.Template entries (single-chain mmCIF plus a
 * query->template residue map) and writes them into each protein chain of the
 * MSA JSON, so the colabfold2 backends predict with templates when they are
 * available.
 */
process AF2_TEMPLATES_TO_AF3_JSON {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/jscgh/proteinfold_colabfold:af3-preview-e9e8bd0"

    input:
    tuple val(meta), path(msa_json), path(raw_dir)
    path(template_mmcif_dir, stageAs: 'template_mmcif_dir')
    path(obsolete_dat, stageAs: 'obsolete.dat')
    val max_template_date

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    af2_templates_to_af3_json.py \\
        ${msa_json} \\
        --raw_dir raw_dir \\
        --template_mmcif_dir template_mmcif_dir \\
        --obsolete obsolete.dat \\
        --output ${meta.id}.json \\
        --max-template-date ${max_template_date} \\
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
