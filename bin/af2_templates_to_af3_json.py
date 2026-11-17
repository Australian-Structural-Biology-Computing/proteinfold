#!/usr/bin/env python3
"""Populate the templates of an AlphaFold3-format MSA JSON from alphafold2 hits.

RUN_ALPHAFOLD2_MSA (the alphafold2 multimer data pipeline) searches templates
per chain with hmmsearch against pdb_seqres and writes <chain>/pdb_hits.sto
beside the MSA outputs, plus alphafold's chain_id_map.json. This script turns
those hits into alphafold3 folding_input.Template entries - a single-chain
mmCIF and a query->template residue map - reusing colabfold's build_templates,
the same code path the colabfold AF3 backend applies to its own template
results. Each protein chain of the JSON gets:

  "templates": [ {"mmcif": ..., "queryIndices": [...], "templateIndices": [...]}, ... ]

Hits are filtered like alphafold's template featurizer: obsolete PDB entries
are remapped through obsolete.dat (dropped when fully removed), the align
ratio must exceed 0.1, templates shorter than 10 residues and near-duplicate
subsequences of the query are dropped, structures released after
--max-template-date are dropped, and at most --max-templates hits per chain
are kept (best first, in hmmsearch order). Chains without usable hits keep
"templates": [].
"""
import argparse
import datetime
import json
import os

MIN_ALIGN_RATIO = 0.1
MIN_TEMPLATE_LENGTH = 10
MAX_SUBSEQUENCE_RATIO = 0.95


def read_text(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def discover_msa_dir(path):
    """Locate the MSA output directory under a staged dir.

    The multimer data pipeline writes a chain_id_map.json; the monomer
    pipeline writes the alignments at the MSA directory root instead.
    """
    map_candidates = []
    root_candidates = []
    for root, dirs, files in os.walk(path):
        if root.count(os.sep) - path.rstrip(os.sep).count(os.sep) > 3:
            dirs[:] = []
            continue
        if "chain_id_map.json" in files:
            map_candidates.append(root)
        elif "uniref90_hits.sto" in files or "uniref90_hits.sto.zst" in files:
            root_candidates.append(root)
    if len(map_candidates) == 1:
        return map_candidates[0]
    if not map_candidates and len(root_candidates) == 1:
        return root_candidates[0]
    raise ValueError(
        f"Expected exactly one MSA output directory under {path} "
        f"(chain_id_map.json or uniref90_hits.sto), found "
        f"chain_id_map.json: {map_candidates}, uniref90: {root_candidates}"
    )


def chain_msa_dirs(msa_dir):
    """Map each chain id to its MSA directory (identical chains share one).

    Returns (chain_id_map or None, chain_ids, dirs); chain_id_map is None for
    the monomer data pipeline layout.
    """
    chain_map_path = os.path.join(msa_dir, "chain_id_map.json")
    if not os.path.isfile(chain_map_path):
        # Monomer data pipeline: single chain, alignments at the MSA dir root.
        return None, ["A"], {"A": msa_dir}
    with open(chain_map_path, encoding="utf-8") as handle:
        chain_id_map = json.load(handle)
    chain_ids = sorted(chain_id_map)
    first_by_sequence = {}
    dirs = {}
    for chain_id in chain_ids:
        sequence = chain_id_map[chain_id]["sequence"]
        first_by_sequence.setdefault(sequence, chain_id)
        dirs[chain_id] = os.path.join(
            msa_dir, first_by_sequence[sequence]
        )
    return chain_id_map, chain_ids, dirs


def parse_obsolete(path):
    """Port of alphafold.data.templates._parse_obsolete."""
    obsolete = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("OBSLTE"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            # OBSLTE date replaced_by ... replaced pdb_id
            obsolete[fields[-1]] = fields[2]
    return obsolete


def template_hits(path, query_sequence):
    """Parse a template search output into alphafold TemplateHit objects.

    The multimer data pipeline writes pdb_hits.sto (hmmsearch against
    pdb_seqres); the monomer pipeline writes pdb_hits.hhr (hhsearch against
    pdb70).
    """
    from alphafold.data import parsers

    if path.endswith(".sto"):
        a3m_string = parsers.convert_stockholm_to_a3m(
            read_text(path), remove_first_row_gaps=False
        )
        return parsers.parse_hmmsearch_a3m(
            query_sequence=query_sequence, a3m_string=a3m_string,
            skip_first=False,
        )
    return parsers.parse_hhr(read_text(path))


def passes_prefilter(hit, query_sequence):
    """Port of alphafold's _assess_hhsearch_hit checks that need no mmCIF."""
    aligned_cols = hit.aligned_cols
    align_ratio = aligned_cols / len(query_sequence)
    if align_ratio <= MIN_ALIGN_RATIO:
        return "align ratio %.2f" % align_ratio
    template_sequence = hit.hit_sequence.replace("-", "")
    length_ratio = len(template_sequence) / len(query_sequence)
    if (
        template_sequence in query_sequence
        and length_ratio > MAX_SUBSEQUENCE_RATIO
    ):
        return "duplicate subsequence of query"
    if len(template_sequence) < MIN_TEMPLATE_LENGTH:
        return "template too short (%d)" % len(template_sequence)
    return None


def release_date_after_cutoff(cif_path, cutoff):
    """True when the structure was released after the cutoff date."""
    from alphafold.data import mmcif_parsing

    pdb_id = os.path.basename(cif_path).split(".")[0]
    parsing_result = mmcif_parsing.parse(
        file_id=pdb_id, mmcif_string=read_text(cif_path)
    )
    if parsing_result.mmcif_object is None:
        return "unparsable mmCIF"
    release_date = parsing_result.mmcif_object.header.get("release_date")
    if not release_date:
        return None  # no date known; keep the template
    try:
        released = datetime.datetime.strptime(release_date, "%Y-%m-%d")
    except (ValueError, TypeError):
        return None
    if released > cutoff:
        return "released %s after cutoff %s" % (release_date, cutoff)
    return None


def select_hits(hits, query_sequence, obsolete, mmcif_dir, cutoff, max_templates):
    """Pick the best valid hits: obsolete/date/prefilter filtered, hmmsearch order."""
    selected = []
    seen = set()
    for hit in hits:
        if len(selected) >= max_templates:
            break
        pdb_id, _, chain_id = hit.name.partition("_")
        chain_id = (chain_id or "A").split()[0]
        key = (pdb_id.lower(), chain_id)
        if key in seen:
            continue
        reason = passes_prefilter(hit, query_sequence)
        if reason:
            continue
        if pdb_id in obsolete:
            if obsolete[pdb_id] is None:
                continue
            pdb_id = obsolete[pdb_id]
            key = (pdb_id.lower(), chain_id)
            if key in seen:
                continue
        cif_path = os.path.join(mmcif_dir, f"{pdb_id.lower()}.cif")
        if not os.path.isfile(cif_path):
            continue
        reason = release_date_after_cutoff(cif_path, cutoff)
        if reason:
            continue
        seen.add(key)
        selected.append(hit)
    return selected


def templates_for_chain(chain_msa_dir, sequence, obsolete, mmcif_dir, cutoff,
                        max_templates):
    """Build the alphafold3 template entries for one chain."""
    from colabfold.alphafold3.templates import build_templates

    hits = None
    for name in ("pdb_hits.sto", "pdb_hits.hhr"):
        candidate = os.path.join(chain_msa_dir, name)
        if os.path.isfile(candidate):
            hits = template_hits(candidate, sequence)
            break
    if hits is None:
        return []
    selected = select_hits(
        hits, sequence, obsolete, mmcif_dir, cutoff, max_templates
    )
    if not selected:
        return []
    templates = build_templates(
        selected, sequence, mmcif_dir, max_templates=max_templates
    )
    entries = []
    for template in templates:
        mapping = sorted(template.query_to_template_map.items())
        entries.append({
            "mmcif": template.mmcif,
            "queryIndices": [query for query, _ in mapping],
            "templateIndices": [target for _, target in mapping],
        })
    return entries


def main():
    parser = argparse.ArgumentParser(
        description="Add alphafold2 template hits to an AlphaFold3-format MSA JSON"
    )
    parser.add_argument("json_in")
    parser.add_argument("--raw_dir", required=True,
                        help="RUN_ALPHAFOLD2_MSA raw output bundle")
    parser.add_argument("--template_mmcif_dir", required=True,
                        help="Directory of template mmCIF files (<pdb_id>.cif)")
    parser.add_argument("--obsolete", required=True,
                        help="PDB obsolete.dat mapping")
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-template-date", required=True,
                        help="Latest template release date to keep")
    parser.add_argument("--max-templates", type=int, default=20)
    args = parser.parse_args()

    cutoff = datetime.datetime.strptime(args.max_template_date, "%Y-%m-%d")
    obsolete = parse_obsolete(args.obsolete)
    msa_dir = discover_msa_dir(args.raw_dir)
    chain_id_map, chain_ids, chain_dirs = chain_msa_dirs(msa_dir)

    document = json.loads(read_text(args.json_in))
    proteins = [
        entity["protein"]
        for entity in document.get("sequences", [])
        if isinstance(entity, dict) and isinstance(entity.get("protein"), dict)
    ]
    if len(proteins) != len(chain_ids):
        raise ValueError(
            f"JSON has {len(proteins)} protein chain(s) but chain_id_map.json "
            f"has {len(chain_ids)}: {chain_ids}"
        )

    for protein, chain_id in zip(proteins, chain_ids):
        sequence = protein["sequence"]
        if chain_id_map is not None and chain_id_map[chain_id]["sequence"] != sequence:
            raise ValueError(
                f"Chain {chain_id} sequence mismatch between the JSON and "
                "chain_id_map.json"
            )
        protein["templates"] = templates_for_chain(
            chain_dirs[chain_id], sequence, obsolete, args.template_mmcif_dir,
            cutoff, args.max_templates,
        )

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(document, indent=2) + "\n")


if __name__ == "__main__":
    main()
