process MULTIFASTA_TO_COLABFOLD_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: query_fasta
    tuple val("${task.process}"), val('sed'), eval('sed --version 2>&1 | sed -n "s/^.*GNU sed) //p"'), emit: versions_sed, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Colon-join every record of the sample FASTA into a single query so
    # complexes are searched (and paired) as one entry, matching the
    # AlphaFold3 JSON extraction. colabfold_search treats a sequence
    # containing ':' as a complex.
    printf '>${meta.id}\\\\n' > ${meta.id}.fasta
    awk '/^>/ { if (records) printf ":"; records++; next } { printf "%s", \\\$0 } END { printf "\\\\n" }' ${fasta} >> ${meta.id}.fasta
    """

    stub:
    """
    printf '>${meta.id}\\\\nstub\\\\n' > ${meta.id}.fasta
    """
}
