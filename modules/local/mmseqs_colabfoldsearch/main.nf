process MMSEQS_COLABFOLDSEARCH_BATCH {
    tag "search_${query_fasta.simpleName}"
    label 'process_high_memory'
    label 'process_high'
    container "ghcr.io/tlitfin/proteinfold_colabfold:af3-preview-38683c2"

    input:
    tuple val(samples), path(query_fasta) // one concatenated ColabFold query (msa_batch_size records)
    path ('db/*')
    path ('uniref30/*')

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

    """

    stub:
    // One dummy a3m/JSON pair per record of the concatenated query; the
    // sample ids come from the record headers, sanitised the way
    // colabfold_search names its results.
    """
    mkdir results
    python3 - ${query_fasta} <<'PY'
import json, os, re, sys

records, header, seq = [], None, []
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith('>'):
        if header is not None:
            records.append((header, ''.join(seq)))
        header, seq = line[1:], []
    elif line:
        seq.append(line)
if header is not None:
    records.append((header, ''.join(seq)))

for raw_id, sequence in records:
    sample_id = re.sub(r'[^A-Za-z0-9_.-]', '_', raw_id)
    with open(f'results/{sample_id}.a3m', 'w') as handle:
        handle.write('>stub\\n')
    document = {
        'dialect': 'alphafold3',
        'version': 2,
        'name': sample_id,
        'modelSeeds': [1],
        'sequences': [
            {
                'protein': {
                    'id': ['A'],
                    'sequence': sequence or 'ACDE',
                    'unpairedMsa': '>' + (raw_id or 'query') + '\\n' + (sequence or 'ACDE') + '\\n',
                }
            }
        ],
    }
    with open(f'results/{sample_id}.json', 'w') as handle:
        json.dump(document, handle)
PY
    """
}
