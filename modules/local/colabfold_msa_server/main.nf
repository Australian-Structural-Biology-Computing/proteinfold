process COLABFOLD_MSA_SERVER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/jscgh/proteinfold_mmseqs_colabfoldsearch:2.1.0dev"

    input:
    tuple val(meta), path(input_csv)

    output:
    tuple val(meta), path("results/*.a3m"), emit: a3m
    tuple val("${task.process}"), val('colabfold'), eval("pip list | grep -i \"^colabfold \" | awk '{print \\\$2}' 2>/dev/null || echo \"unknown\""), emit: versions_colabfold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python ${moduleDir}/run_msa_server.py $args ${input_csv} results/
    test "\$(find results -maxdepth 1 -type f -name '*.a3m' | wc -l)" -eq 1 || { echo "Expected one ColabFold MSA server output per sample" >&2; exit 1; }
    for msa in results/*.a3m; do test -s "\$msa" || { echo "Empty MSA output: \$msa" >&2; exit 1; }; done
    """

    stub:
    """
    mkdir results
    printf '#4\t1\n>101\nACDE\n' > results/${meta.id}.a3m
    """
}
