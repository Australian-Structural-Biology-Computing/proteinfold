include { batchTag } from '../../../subworkflows/local/utils_nfcore_proteinfold_pipeline'

process MMSEQS_COLABFOLDSEARCH_BATCH {
    tag "${batchTag(batch_meta)}"
    label 'process_high_memory'
    label 'process_high'
    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    input:
    tuple val(batch_meta), path(fastas)
    path ('db/*')
    path ('uniref30/*')
    val cache_dir

    output:
    tuple val(batch_meta), path("results/*.a3m"), emit: a3m
    tuple val(batch_meta), path("results/*.json"), emit: json
    tuple val("${task.process}"), val('colabfold_search'), eval("pip list | grep \"^colabfold\" | head -1 | awk '{print \\\$2}' 2>/dev/null || echo \"unknown\""), emit: versions_colabfold_search, topic: versions
    tuple val("${task.process}"), val('mmseqs'), eval("mmseqs version 2>/dev/null | head -1"), emit: versions_mmseqs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local MMSEQS_COLABFOLDSEARCH_BATCH module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def dir = cache_dir ?: ''
    def cache_id = task.ext.msa_cache_id ?: ''
    def tool_id = task.ext.msa_cache_tool_id ?: ''
    def quotedDir = "'${dir.toString().replace("'", "'\\''")}'"
    def quotedId = "'${cache_id.toString().replace("'", "'\\''")}'"
    def quotedToolId = "'${tool_id.toString().replace("'", "'\\''")}'"
    def quotedArgs = "'${args.toString().replace("'", "'\\''")}'"
    """
    for f in uniref30/*; do
        [ -e "db/\$(basename \$f)" ] || ln -sf \$(realpath \$f) db/\$(basename \$f)
    done
    # Every input is a colon-delimited FASTA holding exactly one record whose
    # header is the sample id (multifasta_to_colabfold_fasta or the AF3 JSON
    # extractor), so batching is a plain concatenation and colabfold_search
    # names each result after the record header.
    cat ${fastas} > query.fasta
    colabfold_search $args --threads $task.cpus query.fasta ./db --af3-json results/
    for f in ${fastas}; do
        id=\$(basename "\$f")
        id="\${id%.*}"
        test -s "results/\${id}.a3m" || { echo "Missing A3M output for \$id" >&2; exit 1; }
        test -s "results/\${id}.json" || { echo "Missing JSON output for \$id" >&2; exit 1; }
        if [ -n ${quotedDir} ]; then
            msa_cache.py store "\$f" ${quotedId} ${quotedToolId} ${quotedArgs} ${quotedDir} "results/\${id}.json" "results/\${id}.a3m"
        fi
    done
    """

    stub:
    """
    mkdir results
    for f in ${fastas}; do
        id=\$(basename "\$f"); id="\${id%.*}"
        printf '>stub\n' > results/\${id}.a3m
        python3 - "\$f" "\$id" <<'PY'
import json, sys
path, sample_id = sys.argv[1], sys.argv[2]
records, header, seq = [], None, []
for line in open(path):
    line = line.strip()
    if line.startswith('>'):
        if header is not None:
            records.append((header, ''.join(seq)))
        header, seq = line[1:], []
    elif line:
        seq.append(line)
if header is not None:
    records.append((header, ''.join(seq)))
document = {
    'dialect': 'alphafold3',
    'version': 2,
    'name': sample_id,
    'modelSeeds': [1],
    'sequences': [
        {
            'protein': {
                'id': [chr(ord('A') + i) for i in range(len(records))] or ['A'],
                'sequence': ':'.join(s for _, s in records) or 'ACDE',
                'unpairedMsa': '>' + (records[0][0] if records else 'query') + '\\n' + (':'.join(s for _, s in records) or 'ACDE') + '\\n',
            }
        }
    ],
}
with open(f'results/{sample_id}.json', 'w') as handle:
    json.dump(document, handle)
PY
    done
    """
}
