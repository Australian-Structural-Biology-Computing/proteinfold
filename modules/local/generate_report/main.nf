include { batchTag } from '../../../subworkflows/local/utils_nfcore_proteinfold_pipeline'

process GENERATE_REPORT_BATCH {
    tag   "${batchTag(batch_meta)}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/241f0746484727a3633f544c3747bfb77932e1c8c252e769640bd163232d9112/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pip_plotly:35975fa0fc54b2d3' }"

    input:
    tuple val(batch_meta), path(structures), path(msa_files), path(metrics_stores), val(manifest), val(struct_ext)
    path(template)

    output:
    tuple val(batch_meta), path ("*report.html")     , emit: report
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('generate_report.py'), eval("python3 --version | sed 's/Python //g'"), emit: versions_generator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model = (batch_meta instanceof Collection ? batch_meta[0] : batch_meta).with { m -> m.mode_short ?: m.model }
    // Template written at column 0: the multi-line ${manifest} insertion
    // defeats Nextflow's script indentation stripping, which previously left
    // heredoc terminators indented (breaking bash) in the batch inference
    // script. See the matching fix in colabfold2_inference.
    """
cat > sample_manifest.tsv <<BATCH_MANIFEST
${manifest}
BATCH_MANIFEST

export PYTHONDONTWRITEBYTECODE=1

generate_report.py \
    --manifest sample_manifest.tsv \
    --struct-ext ${struct_ext} \
    --type ${model} \
    --html_template ${template} \
    --output_dir . \
    ${args}
"""

    stub:
    def model = (batch_meta instanceof Collection ? batch_meta[0] : batch_meta).with { m -> m.mode_short ?: m.model }
    def stub_files = (batch_meta instanceof Collection ? batch_meta.collect() : [batch_meta]).collect { meta -> """
    touch ${meta.id}_${meta.mode_short}_${meta.msa_short}_report.html
    touch ${meta.id}_${meta.mode_short}_${meta.msa_short}_seq_coverage.png
    touch ${meta.id}_coverage_LDDT.html
    """ }.join('\n')
    """
${stub_files}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_report.py: \$(python3 --version)
    END_VERSIONS
    """
}
