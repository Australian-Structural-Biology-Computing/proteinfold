#!/usr/bin/env python3

import os
import argparse
import json
import sys
import warnings
#import torch moved to a conditional import since too bulky import if not used
import numpy as np
import csv
import re
from utils import plddt_from_struct_b_factor, get_chain_ids

# TODO: Issue #309, make into a proper separate process, it its own module so that dependencies can be managed better
# TODO: Need a sense of ranking, so that metrics can be traced back to correct model structure, even if they're not in sequential order. The enumerates() here are not sufficient.
#       Needs to be program-dependent, (see item below).
# TODO: look into have a --prog argument that could set filenames etc, logically seperate it?
# {name}_{prog}_{metric}.tsv might be easier for MultiQC to parse a complex workdir, than without the .prog
# TODO: read --prog from ${meta.model} in the NextFlow pipes. This also allows case switching in a proper EXTRACT_METRICS process.
# E.g. in main.nf of EXTACT_METRICS process, we could have:
# match ${meta.mode}:
#     case 'alphafold2':
#        ...
#        ...
#...
# ^ overwrought with duplication, but can catch program specific weirdness, and lower barrier to adding new programs in the future.

# TODO: Chain-wise iPTM since the relevant interface might not always be the average of all.
# Would complete Issue #308
# Proposed format is pair-interfaces in rows, structure inference number in cols: https://github.com/nf-core/proteinfold/pull/312#issuecomment-2917709432
# KR - changed to have both sides of the matrix, because it's not symmetrical (see comment in Issue #306)

# Mapping of characters to integers for MSA parsing.
# 20 is for unknown characters, and 21 is for gaps.
AA_to_int = {
    "A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9,
    "M": 10, "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19,
   ".": 20, "-": 21
}

def a3m_to_int(a3m_file):
    """
    Convert an A3M MSA representation into an integer representation (0-21).
    """
    with open(a3m_file, "r") as f:
        msa = f.read()

    # Convert each sequence in the MSA
    int_sequences = []
    for idx, line in enumerate(msa.splitlines()):
        if idx == 0 and not line.startswith(">"):  # If there's an additional header (non-FASTA) skip it. E.g ColabFold
            continue

        if not line.startswith(">"):  # Ignore header lines
            filtered_line = ''.join(char for char in line if not char.islower()) # Remove inserts (lower-case chars) in a3m
            int_sequence = [AA_to_int.get(char.upper(), 20) for char in filtered_line]
            int_sequences.append(int_sequence)

    int_sequences_array = np.array(int_sequences, dtype=object)
    return int_sequences_array

def format_msa_rows(msa_data):
    return [[str(x) for x in val] for val in msa_data]

def format_pae_rows(pae_data):
    return [[f"{num:.4f}" for num in row] for row in pae_data]

def format_pair_score_rows(pair_score_entries, pair_labels=None):
    if pair_labels is None:
        pair_labels = sorted({label for score_values in pair_score_entries.values() for label, _ in score_values})

    rows = [[""] + pair_labels]
    for model_idx, score_values in sorted(pair_score_entries.items(), key=lambda item: sort_model_label(item[0])):
        score_map = {label: value for label, value in score_values}
        rows.append([model_idx] + [f"{score_map[label]:.4f}" if label in score_map else "n/a" for label in pair_labels])

    return [list(row) for row in zip(*rows)]

def write_tsv(file_path, rows):
    with open(file_path, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerows(rows)

def sort_model_label(label):
    try:
        return (0, int(label))
    except (TypeError, ValueError):
        return (1, str(label))

def infer_model_rank(file_path):
    normalized_path = file_path.replace(os.sep, "/")
    rank_patterns = [
        r"ranked_(\d+)",
        r"_rank_(\d+)",
        r"-rank(\d+)(?:/|$)",
        r"_model_(\d+)",
    ]

    for pattern in rank_patterns:
        match = re.search(pattern, normalized_path)
        if match:
            return int(match.group(1))

    return None

def sort_paths_by_rank(paths):
    def sort_key(path):
        rank = infer_model_rank(path)
        if rank is None:
            warnings.warn(f"Unable to infer model rank from path: {path}; falling back to basename sort")
            return (1, os.path.basename(path))
        return (0, rank, os.path.basename(path))

    return sorted(paths, key=sort_key)

def colabfold_chainwise_rows(data):
    """Read the chainwise metrics surfaced in a colabfold scores JSON.

    The chain-pair metrics come from the scores JSON: the AF3-family
    backends write their confidence head's chain_pair_iptm values and
    per-chain pTM (via the colabfold_batch_af3_chainwise.py wrapper), the
    AlphaFold2 backends write the --calc-extra-ptm pairwise metrics. Both
    report one entry per unordered pair: the on-diagonal (A:A) rows are
    the per-chain pTM and the off-diagonal (A:B) rows the pairwise
    interface ipTM. Returns an empty list for monomer models, which have
    no interface pairs.
    """
    entries = []
    per_chain_ptm = data.get("per_chain_ptm") or {}
    for chain_id, value in per_chain_ptm.items():
        try:
            entries.append((f"{chain_id}:{chain_id}", float(value)))
        except (TypeError, ValueError):
            continue
    for pair_label, value in (data.get("pairwise_iptm") or {}).items():
        chains = pair_label.split("-")
        if len(chains) != 2:
            continue
        try:
            entries.append((f"{chains[0]}:{chains[1]}", float(value)))
        except (TypeError, ValueError):
            continue
    # No off-diagonal pairs means a monomer: nothing beyond the global
    # scores is worth emitting.
    has_pairs = any(label.split(":")[0] != label.split(":")[1] for label, _ in entries)
    return entries if has_pairs else []

def _scores_json_is_structure_only(scores_files):
    if not scores_files:
        return False
    parsed = []
    for fn in scores_files:
        try:
            with open(fn) as f:
                parsed.append((fn, json.load(f)))
        except (OSError, json.JSONDecodeError):
            parsed.append((fn, None))
    if all(data is None for _, data in parsed):
        raise SystemExit(
            f"None of the {len(scores_files)} paired scores JSON file(s) could be read "
            f"(first failure: {parsed[0][0]}). The extraction inputs are broken, so the "
            "run cannot be classified as structure-only or confidence; refusing to "
            "emit a potentially misleading pLDDT table. Check that the scores JSONs "
            "exist, are readable and contain valid JSON."
        )
    for _, data in parsed:
        if data is None:
            return False
        if "ptm" in data:
            return False
        ranking_score = data.get("ranking_score")
        if ranking_score is not None:
            try:
                if not np.isnan(float(ranking_score)):
                    return False
            except (TypeError, ValueError):
                return False
    return True

def _model_type_from_scores_path(scores_path):
    match = re.search(r"_scores_rank_\d+_(.+)_seed_\d+_sample_\d+\.json$", os.path.basename(scores_path))
    return match.group(1) if match and match.group(1) else "unknown"

def extract_structs_plddt_to_tsv(name, structures, scores_files=None):
    """
    Write out a tsv file contain pLDDTs for reading by MultiQC in nf-core/proteinfold
    Uses utils function with BioPython PDB package to extract residue pLDDT values from the b-factor column.
    """
    sorted_structures = sort_paths_by_rank(structures)
    plddt_cols = [plddt_from_struct_b_factor(structure) for structure in sorted_structures]
    res_counts = [len(plddt_col) for plddt_col in plddt_cols]

    if len(set(res_counts)) != 1:
        raise ValueError("Not all structures have the same number of residues!")

    structure_only = all(np.all(plddt_col <= 0) for plddt_col in plddt_cols) and _scores_json_is_structure_only(scores_files)
    if structure_only:
        model_types = sorted({_model_type_from_scores_path(fn) for fn in scores_files})
        print(f"{name}: structure-only model ({', '.join(model_types)}): pLDDT/PAE unavailable (no confidence head)")

    rank_names = []
    for idx, structure in enumerate(sorted_structures):
        rank = infer_model_rank(structure)
        if rank is None:
            warnings.warn(f"Unable to infer model rank from path: {structure}; falling back to index {idx}")
        rank_names.append(f"rank_{rank}" if rank is not None else f"rank_{idx}")
    # Create header as the first row
    plddt_rows =  [["Positions"] + rank_names]
    if structure_only:
        plddt_rows.extend([[res_id] + ["n/a"] * len(plddt_cols) for res_id in range(res_counts[0])])
    else:
        res_id_col = list(range(len(plddt_cols[0])))
        plddt_rows.extend([list(row) for row in zip(res_id_col, *plddt_cols)])  # Combine lists column-wise to make rows
    write_tsv(f"{name}_plddt.tsv", plddt_rows)

def read_paired_a3m(name, a3m_file):
    msa_rows = a3m_to_int(a3m_file)
    write_tsv(f"{name}_msa.tsv", format_msa_rows(msa_rows))

def read_colabfold_metrics(name, colabfold_metrics_files):
    ptm_rows = []
    iptm_rows = []
    ipsae_rows = []
    chainwise = {}
    chainwise_ipsae = {}
    for fn in colabfold_metrics_files:
        with open(fn) as f:
            data = json.load(f)
        rank_id = infer_model_rank(fn)
        if rank_id is None:
            raise ValueError(f"Unable to infer ColabFold rank from metrics filename: {fn}")
        if "pae" in data:
            write_tsv(f"{name}_{rank_id}_pae.tsv", format_pae_rows(data["pae"]))
        if "ptm" in data:
            ptm_rows.append((f"{rank_id}", data["ptm"]))
        if "iptm" in data:
            iptm_rows.append((f"{rank_id}", data["iptm"]))
        # Global ipSAE: colabfold's native per-chain-pair ipSAE dict from
        # get_interface_scores (the AlphaFold2 backends run it natively;
        # the AF3-family backends get it via
        # colabfold_batch_af3_chainwise.py).
        interface_ipsae = data.get("ipsae")
        if isinstance(interface_ipsae, dict) and interface_ipsae:
            ipsae_rows.append((f"{rank_id}", f"{max(interface_ipsae.values()):.3f}"))
            # The native per-chain-pair ipSAE dict is asymmetric (the
            # alignment runs from each chain against the other), so both
            # orientations are kept in the full matrix.
            chainwise_ipsae[rank_id] = [
                (f"{label.replace('-', ':')}", value) for label, value in interface_ipsae.items()
            ]
        chainwise_entries = colabfold_chainwise_rows(data)
        if chainwise_entries:
            chainwise[rank_id] = chainwise_entries
    if len(ptm_rows)>0:
        write_tsv(f"{name}_ptm.tsv", sorted(ptm_rows, key = lambda x: x[0]))
    if len(iptm_rows)>0:
        write_tsv(f"{name}_iptm.tsv", sorted(iptm_rows, key = lambda x: x[0]))
    if len(ipsae_rows)>0:
        write_tsv(f"{name}_ipsae.tsv", sorted(ipsae_rows, key = lambda x: x[0]))
    if len(chainwise)>0:
        write_tsv(f"{name}_chainwise.tsv", format_pair_score_rows(dict(sorted(chainwise.items(), key=lambda x: x[0]))))
    if len(chainwise_ipsae)>0:
        write_tsv(f"{name}_chainwise_ipsae.tsv", format_pair_score_rows(dict(sorted(chainwise_ipsae.items(), key=lambda x: x[0]))))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paired_a3m", dest="paired_a3m", required=False) # For reading the ColabFold MSA format
    parser.add_argument("--colabfold_metrics_files", required=False, nargs="+")
    parser.add_argument("--structs", dest="structs", required=False, nargs="+")
    parser.add_argument("--name", default="untitled", dest="name") # might need a --name $meta.id
    args = parser.parse_args()

    if args.paired_a3m:
        read_paired_a3m(args.name, args.paired_a3m)
    if args.structs:
        extract_structs_plddt_to_tsv(args.name, args.structs, args.colabfold_metrics_files)
    if args.colabfold_metrics_files:
        read_colabfold_metrics(args.name, args.colabfold_metrics_files)

if __name__ == "__main__":
    main()
