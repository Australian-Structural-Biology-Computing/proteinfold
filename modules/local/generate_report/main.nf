include { batchTag } from '../../../subworkflows/local/utils_nfcore_proteinfold_pipeline'

process GENERATE_REPORT_BATCH {
    tag   "${batchTag(batch_meta)}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/241f0746484727a3633f544c3747bfb77932e1c8c252e769640bd163232d9112/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pip_plotly:35975fa0fc54b2d3' }"

    input:
    tuple val(batch_meta), path(structures), path(msa_files), path(pae_files), path(iptm_files), path(ipsae_files), path(chainwise_files), path(chainwise_ipsae_files)
    path(template)

    output:
    tuple val(batch_meta), path("*report.html"), emit: report
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('generate_report.py'), eval("python3 --version | sed 's/Python //g'"), emit: versions_generate_report, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def head = batch_meta instanceof Collection ? batch_meta[0] : batch_meta
    def model = head.mode_short ?: head.model
    // The task files carry the sanitised ids colabfold_batch bakes into its
    // outputs; the report name keeps the raw sample id.
    def raw_ids = (batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta])
        .collect { meta -> safe(meta.id) }.join(' ')
    def sample_ids = (batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta])
        .collect { meta -> (meta.id as String).replaceAll(/[^A-Za-z0-9_.-]/, '_') }
    if (sample_ids.toSet().size() != sample_ids.size()) {
        error("GENERATE_REPORT_BATCH received sample IDs that collide after sanitisation: ${sample_ids}")
    }
    def struct_ext = (structures instanceof List ? structures[0] : structures).name.endsWith('.pdb') ? 'pdb' : 'cif'
    def ids_content = sample_ids.collect { safe(it) }.join(' ')

    """
printf '%s\\n' ${ids_content} > sample_ids.txt
printf '%s\\n' ${raw_ids} > raw_ids.txt

export PYTHONDONTWRITEBYTECODE=1

while IFS= read -r id && IFS= read -r raw_id <&3; do
    generate_report.py \\
        --type ${model} \\
        --tag ${model}_${head.msa_short ?: 'msa'} \\
        --msa "\${id}_colabfold2_msa.tsv" \\
        --pae "\${id}_1_pae.tsv" \\
        --iptm "\${id}_iptm.tsv" \\
        --ipsae "\${id}_ipsae.tsv" \\
        --chainwise "\${id}_chainwise.tsv" \\
        --chainwise_ipsae "\${id}_chainwise_ipsae.tsv" \\
        --pdb "\${id}"_rank_*.${struct_ext} \\
        --html_template ${template} \\
        --output_dir ./ \\
        --name "\${raw_id}" \\
        $args || exit 1
done < sample_ids.txt 3< raw_ids.txt
"""

    stub:
    def head = batch_meta instanceof Collection ? batch_meta[0] : batch_meta
    def stub_files = (batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]).collect { meta ->
        def tag = "${meta.mode_short ?: meta.model}_${meta.msa_short ?: 'msa'}"
"    touch ${meta.id}_${tag}_report.html\n"
    }.join()
    """
${stub_files}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_report.py: \$(python3 --version)
    END_VERSIONS
    """
}

def safe(name) {
    "'" + (name as String).replaceAll("'", "'\\''") + "'"
}
