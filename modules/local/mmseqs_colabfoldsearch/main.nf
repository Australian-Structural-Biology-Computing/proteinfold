process MMSEQS_COLABFOLDSEARCH_BATCH {
    tag "search_${query_fasta.simpleName}"
    label 'process_high_memory'
    label 'process_high'
    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    input:
    tuple val(samples), path(query_fasta) // one concatenated ColabFold query (msa_batch_size records)
    path ('db/*')
    path ('uniref30/*')
    val cache_dir

    output:
    tuple val(samples), path("results/*.a3m"), emit: a3m
    tuple val(samples), path("results/*.json"), emit: json
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
    colabfold_search $args --threads $task.cpus ${query_fasta} ./db --af3-json results/

    # colabfold_search names each result after the record header, sanitised
    # to [A-Za-z0-9_.-]; the sample ids come from the query's records, not
    # from the task metadata.
    grep '^>' ${query_fasta} | sed 's/^>//' | sed 's/[^A-Za-z0-9_.-]/_/g' | while read -r id; do
        test -s "results/\${id}.a3m" || { echo "Missing A3M output for \$id" >&2; exit 1; }
        test -s "results/\${id}.json" || { echo "Missing JSON output for \$id" >&2; exit 1; }
    done

    if [ -n ${quotedDir} ]; then
        # The store keys on the per-sample query content, so split the
        # concatenated query back into transient per-record files; they are
        # removed once stored.
        mkdir keys
        awk '/^>/ { if (f) close(f); h = substr(\$0, 2); n = h; gsub(/[^A-Za-z0-9_.-]/, "_", n); f = "keys/" n ".fasta" } f { print > f }' ${query_fasta}
        for f in keys/*.fasta; do
            id=\$(basename "\$f")
            id="\${id%.*}"
            msa_cache.py store "\$f" ${quotedId} ${quotedToolId} ${quotedArgs} ${quotedDir} "results/\${id}.json" "results/\${id}.a3m"
        done
        rm -rf keys
    fi
    """

    stub:
    """
    mkdir results keys
        awk '/^>/ { if (f) close(f); h = substr(\$0, 2); n = h; gsub(/[^A-Za-z0-9_.-]/, "_", n); f = "keys/" n ".fasta" } f { print > f }' ${query_fasta}
    for f in keys/*.fasta; do
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
    rm -rf keys
    """
}
