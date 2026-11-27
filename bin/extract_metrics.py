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

def parse_a3m_cardinality_header(first_line):
    """Parse a ColabFold complex-a3m header '#len1,...,lenN\tcard1,...,cardN'.

    ColabFold writes complex a3m files with every MSA row collapsed to one
    span per unique entity (colabfold.input.msa_to_str resets cardinality to
    1 before building the rows); the copy counts live only in this header,
    and ColabFold's own reader expands them. Returns (lengths, cardinalities)
    or None for plain single-chain a3m files and unrecognised headers.
    """
    if not first_line.startswith("#"):
        return None
    fields = first_line[1:].split("\t")
    if len(fields) < 2:
        return None
    try:
        lengths = [int(x) for x in fields[0].split(",")]
        cards = [int(x) for x in fields[1].split(",")]
    except ValueError:
        return None
    if not lengths or len(lengths) != len(cards):
        return None
    if any(length <= 0 for length in lengths) or any(card <= 0 for card in cards):
        return None
    return lengths, cards

def read_structure_polymer_sequences(struct_file):
    """Ordered (chain id, one-letter sequence) for a structure's polymer chains.

    Order matters: prediction outputs arrange chains in the featurised batch
    order, which for AF3-family backends is the sorted-chain order rather than
    the input order, and the MSA coverage plot aligns its columns against this
    order. Waters are skipped and non-standard residues map to X so modified
    polymer chains still match their a3m entity; ligand-only chains (every
    residue a HETATM) are excluded since they carry no MSA columns.
    """
    from Bio import PDB
    if str(struct_file).endswith(".pdb"):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("msa_order", struct_file)
    elif str(struct_file).endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure("msa_order", struct_file)
    else:
        raise ValueError(f"{struct_file} is neither a PDB or mmCIF file!")
    chains = []
    for model in structure:
        for chain in model:
            letters = []
            has_polymer = False
            for residue in chain:
                hetflag = residue.id[0]
                if hetflag == "W":
                    continue
                if hetflag == " ":
                    has_polymer = True
                try:
                    letters.append(PDB.Polypeptide.three_to_one(residue.resname))
                except (ValueError, KeyError):
                    letters.append("X")
            # A chain whose residues are all HETATM records is a ligand (or an
            # exotic single-modified-residue chain, which the sequence match
            # below would reject anyway, falling back to the collapsed MSA).
            if letters and has_polymer:
                chains.append((chain.id, "".join(letters)))
        break  # first model only - predictions carry a single model
    return chains

def order_copies_by_chain_ids(id_spans, structure_chain_ids):
    """Order an MSA's entity copies into a structure's chain order by id.

    The predicted structure arranges chains in the featurised batch order -
    for AF3-family backends the sorted-chain order, which multi-copy inputs
    (entity 'id' lists) interleave relative to the input order the MSA rows
    follow. Each chain id claims the entity copy carrying that id; ligand and
    nucleic chains claim nothing (they have no protein MSA columns). Returns
    the (entity, copy) pairs in structure order, or None when a copy stays
    unmatched or an id repeats - the caller then keeps the collapsed MSA
    unordered, which the viewer safely leaves unlinked rather than
    misaligning bands against the wrong chain.
    """
    claimed = {}
    for chain_id, entity, copy in id_spans:
        if chain_id is None or chain_id in claimed:
            return None
        claimed[chain_id] = (entity, copy)
    order = []
    for chain_id in structure_chain_ids:
        pair = claimed.pop(chain_id, None)
        if pair is not None:
            order.append(pair)
    if claimed or not order:
        return None
    return order


def order_copies_by_sequences(entity_sequences, cards, structure_sequences):
    """Order an MSA's entity copies into a structure's chain order by sequence.

    Used for a3m inputs, which carry no chain ids: each structure polymer
    chain is matched to an unused entity copy by exact one-letter sequence.
    Copies of one entity share a sequence and MSA, so they are
    interchangeable. Returns the (entity, copy) pairs in structure order, or
    None when any chain (or copy) stays unmatched.
    """
    groups = [(entity, copy) for entity, card in enumerate(cards) for copy in range(card)]
    used = set()
    order = []
    for _, chain_sequence in structure_sequences:
        match = next(
            (
                g for g, (entity, _) in enumerate(groups)
                if g not in used and entity_sequences[entity] == chain_sequence
            ),
            None,
        )
        if match is None:
            return None
        used.add(match)
        order.append(groups[match])
    if len(used) != len(groups):
        return None
    return order


class MissingMsaError(ValueError):
    """The job carries no protein MSA at all (e.g. language-model variants)."""


class NonProteinOnlyError(ValueError):
    """The job has no protein entity, so there is nothing to cover."""


def _a3m_text_rows(text):
    """Sequence rows of an embedded a3m text block, in file order."""
    if not text or not text.strip():
        return []
    if "\n" not in text and not text.startswith(">"):
        # A file path rather than embedded text; the merge step always embeds,
        # so treat anything else as absent.
        return []
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            rows.append("")
        elif rows:
            rows[-1] += line
        else:
            rows.append(line)
    return rows


def af3_msa_json_to_int(json_path, structure=None):
    """Build the collapsed MSA coverage matrix from an AlphaFold3-format JSON.

    The JSON is the merged MSA the inference consumed, and - unlike the a3m
    ColabFold writes back out - it names every chain: an entity's 'id' list
    holds one entry per copy. The rows stay collapsed to one span per entity
    (the true merged depth); the copy layout is matched against the predicted
    structure's chain ids and returned as the display span order for the
    report, which splits the copies out after sorting the rows by identity -
    so every copy of a subunit is sorted and coloured from the same
    per-entity row. Proteins only: nucleic entities have no protein MSA
    columns.

    Returns (rows, spans): the integer rows in entity order plus the
    (length, entity) span layout in display column order, or spans=None when
    the chain ids cannot be matched against the structure - the caller then
    writes the plain collapsed matrix, which the report leaves unlinked
    rather than misaligned.
    """
    with open(json_path, encoding="utf-8") as handle:
        document = json.load(handle)
    if isinstance(document, list):
        if len(document) != 1:
            raise ValueError(
                "ColabFold2 requires one AlphaFold3 job per samplesheet row; "
                f"the JSON contains {len(document)} jobs"
            )
        document = document[0]

    entities = []
    for entity in document.get("sequences", []):
        if not isinstance(entity, dict) or "protein" not in entity:
            continue
        protein = entity["protein"]
        sequence = protein.get("sequence") if isinstance(protein, dict) else None
        if not sequence:
            continue
        ids = protein.get("id")
        id_list = list(ids) if isinstance(ids, list) else ([ids] if ids else [])
        entities.append({
            "ids": id_list,
            "sequence": sequence,
            "paired": _a3m_text_rows(protein.get("pairedMsa")),
            "unpaired": _a3m_text_rows(protein.get("unpairedMsa")),
        })
    if not entities:
        raise NonProteinOnlyError(
            f"{json_path}: the job has no protein entity; there is no MSA "
            "coverage to report"
        )
    if not any(entity["paired"] or entity["unpaired"] for entity in entities):
        raise MissingMsaError(
            f"{json_path}: no protein entity carries an MSA"
        )

    lengths = [len(entity["sequence"]) for entity in entities]
    rows = _coverage_rows(entities, lengths, json_path)
    if rows is None:
        raise ValueError(
            f"{json_path}: the MSA rows do not match the declared entity "
            "sequences"
        )
    if not rows:
        raise MissingMsaError(
            f"{json_path}: no protein entity carries an MSA"
        )

    spans = None
    if structure is not None and all(entity["ids"] for entity in entities):
        id_spans = [
            (chain_id, entity, copy)
            for entity, id_list in enumerate(entity["ids"] for entity in entities)
            for copy, chain_id in enumerate(id_list)
        ]
        order = order_copies_by_chain_ids(id_spans, get_chain_ids(structure))
        if order is not None:
            spans = [(lengths[entity], entity) for entity, _ in order]
        else:
            warnings.warn(
                f"{json_path}: could not match every MSA entity copy to a "
                "predicted chain by id; keeping the collapsed per-entity MSA "
                "so the coverage plot stays unlinked instead of misaligned"
            )
    return rows, spans


def _coverage_rows(entities, lengths, json_path):
    """Assemble the collapsed coverage rows, one span per entity.

    Paired blocks are row-aligned across entities; a block that is missing or
    shallower than the deepest one pads as gaps, which is how the featuriser
    reads a chain without a paired MSA. Unpaired rows follow their entity's
    paired block, gaps across the other entities. The report splits entity
    copies out of these rows after sorting, so no copy expansion happens
    here. Returns the rows, [] when no entity carries an MSA at all, or None
    when a row does not match its entity's declared sequence length.
    """

    def stripped(row):
        return "".join(char for char in row if not char.islower())

    paired_depths = [len(entity["paired"]) for entity in entities]
    paired_depth = max(paired_depths) if paired_depths else 0

    rows = []
    for depth in range(paired_depth):
        row = []
        for entity, length in zip(entities, lengths):
            if depth < len(entity["paired"]):
                values = [AA_to_int.get(char.upper(), 20) for char in stripped(entity["paired"][depth])]
            else:
                values = [21] * length
            if len(values) != length:
                return None
            row.extend(values)
        rows.append(row)
    if not rows:
        # No entity carries a paired block: row 0 is still the full query.
        rows.append([
            AA_to_int.get(char.upper(), 20)
            for entity in entities
            for char in entity["sequence"]
        ])
    for entity_index, (entity, length) in enumerate(zip(entities, lengths)):
        for unpaired_row in entity["unpaired"]:
            segment = [AA_to_int.get(char.upper(), 20) for char in stripped(unpaired_row)]
            if len(segment) != length:
                return None
            rows.append(
                [21] * sum(lengths[:entity_index])
                + segment
                + [21] * sum(lengths[entity_index + 1:])
            )
    return rows


def a3m_to_int(a3m_file, structure=None):
    """
    Convert an A3M MSA representation into an integer representation (0-21).

    Complex a3m files written by ColabFold collapse every entity to one span
    and carry the copy counts in a '#lengths\tcardinalities' header. The rows
    stay collapsed (the true merged depth); when a predicted structure is
    available the copy layout is matched against the structure's chains - by
    chain id when the a3m is one the pipeline re-serialised from a named
    input, otherwise by sequence - and returned as the display span order for
    the report, which splits the copies out after sorting the rows by
    identity. Any deviation from the expected layout returns the plain
    collapsed conversion without a layout, which downstream viewers treat as
    unlinked rather than misaligned.

    Returns (rows, spans); spans is None outside the '#'-header case or when
    the copies cannot be matched.
    """
    with open(a3m_file, "r") as f:
        lines = f.read().splitlines()

    # Convert each sequence in the MSA (the '#'-header line, if any, is
    # skipped like any other non-FASTA first line).
    int_sequences = []
    for idx, line in enumerate(lines):
        if idx == 0 and not line.startswith(">"):  # If there's an additional header (non-FASTA) skip it. E.g ColabFold
            continue

        if not line.startswith(">"):  # Ignore header lines
            filtered_line = ''.join(char for char in line if not char.islower()) # Remove inserts (lower-case chars) in a3m
            int_sequence = [AA_to_int.get(char.upper(), 20) for char in filtered_line]
            int_sequences.append(int_sequence)

    int_sequences_array = np.array(int_sequences, dtype=object)

    spans = None
    card = parse_a3m_cardinality_header(lines[0] if lines else "")
    if card is not None and structure is not None and len(int_sequences):
        lengths, cards = card
        if sum(lengths) == len(int_sequences[0]):
            int_to_letter = {v: k for k, v in AA_to_int.items() if k != "."}
            entity_sequences = [
                "".join(int_to_letter.get(value, "X") for value in int_sequences[0][offset:offset + length])
                for offset, length in zip(np.cumsum([0] + lengths[:-1]), lengths)
            ]
            order = order_copies_by_sequences(
                entity_sequences, cards, read_structure_polymer_sequences(structure)
            )
            if order is not None:
                spans = [(lengths[entity], entity) for entity, _ in order]
            else:
                warnings.warn(
                    f"{a3m_file}: could not match every MSA entity copy to a "
                    "predicted chain; keeping the collapsed MSA so the "
                    "coverage plot stays unlinked instead of misaligned"
                )
    return int_sequences_array, spans

def write_msa_tsv(name, rows, spans=None):
    """Write the MSA TSV, with a '#spans' header describing the display layout.

    The rows stay collapsed to one span per entity - the true merged MSA
    depth - and the header carries the per-copy span layout (span lengths and
    their entity indices, in display column order) that the report expands
    after sorting the rows by identity, so every copy of a subunit inherits
    the same sort position and identity colour. The header is only written
    when the copy layout is verified against the predicted structure;
    without it consumers see a plain matrix and the report leaves the
    coverage unlinked.
    """
    with open(f"{name}_msa.tsv", "w") as out_f:
        if spans is not None:
            lengths = ",".join(str(length) for length, _ in spans)
            entities = ",".join(str(entity) for _, entity in spans)
            out_f.write(f"#spans\t{lengths}\t{entities}\n")
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerows(format_msa_rows(rows))

def format_msa_rows(msa_data):
    # A generator: copy-expanded matrices reach tens of millions of cells and
    # materialising them all as strings would need several gigabytes.
    for val in msa_data:
        yield [str(x) for x in val]

def format_pae_rows(pae_data):
    return [[f"{num:.4f}" for num in row] for row in pae_data]

def format_token_chain_rows(token_chain_ids, pae_data):
    """Header rows describing the PAE token space, or [] when unavailable.

    AF3-dialect outputs score one token per polymer residue but one token per
    atom for ligands, so the matrix can cover more positions than the
    structure has residues. The per-token chain ids let the report draw
    correct chain boundaries (ligand chains included) and map tokens back
    onto structure residues. Encoded as run-length pairs so the header stays
    short even for wide matrices; consumers must verify that the run lengths
    sum to the matrix size before trusting it.
    """
    if not isinstance(token_chain_ids, list) or len(token_chain_ids) != len(pae_data):
        return []
    runs = []
    for chain_id in token_chain_ids:
        if runs and runs[-1][0] == chain_id:
            runs[-1][1] += 1
        else:
            runs.append([str(chain_id), 1])
    return [["# token_chains", json.dumps(runs)]]

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

def read_paired_a3m(name, a3m_file, structure=None):
    msa_rows, spans = a3m_to_int(a3m_file, structure=structure)
    write_msa_tsv(name, msa_rows, spans=spans)

def entity_map_from_af3_json(json_path, structure=None):
    """Chain groups per entity from an AlphaFold3-format JSON.

    Every entity kind (protein, ligand, RNA, DNA) carries an 'id' or an id
    list, with one entry per copy of the entity - the grouping the report's
    entity colouring needs, which the predicted structure alone cannot
    provide (it has no notion of which chains are copies). Returns
    [{'label': 'A+C', 'chains': ['A', 'C']}, ...] in input order, each
    entity's chains listed in the predicted structure's chain order so the
    labels read in the order the viewer shows them.
    """
    with open(json_path, encoding="utf-8") as handle:
        document = json.load(handle)
    if isinstance(document, list):
        if len(document) != 1:
            raise ValueError(
                "ColabFold2 requires one AlphaFold3 job per samplesheet row; "
                f"the JSON contains {len(document)} jobs"
            )
        document = document[0]
    entities = []
    for entity in document.get("sequences", []):
        if not isinstance(entity, dict):
            continue
        payload = next((v for v in entity.values() if isinstance(v, dict)), None)
        if payload is None:
            continue
        ids = payload.get("id")
        id_list = list(ids) if isinstance(ids, list) else ([ids] if ids else [])
        if id_list:
            entities.append(id_list)
    if not entities:
        return []
    if structure is not None:
        chain_order = {chain_id: ix for ix, chain_id in enumerate(get_chain_ids(structure))}
        entities = [
            sorted(id_list, key=lambda cid: chain_order.get(cid, len(chain_order)))
            for id_list in entities
        ]
    return [{"label": "+".join(id_list), "chains": list(id_list)} for id_list in entities]


def read_af3_msa_json(name, json_path, structure=None):
    # The entity map is written before the coverage: it only needs the ids,
    # and must survive the no-MSA exits so language-model variants keep their
    # entity colouring too.
    entity_map = entity_map_from_af3_json(json_path, structure=structure)
    with open(f"{name}_entities.json", "w") as handle:
        json.dump({"entities": entity_map}, handle)
    try:
        msa_rows, spans = af3_msa_json_to_int(json_path, structure=structure)
    except MissingMsaError:
        # The module contract: exit 2 leaves an empty MSA TSV behind, which
        # downstream treats the same as a missing MSA.
        sys.exit(2)
    except NonProteinOnlyError:
        sys.exit(3)
    write_msa_tsv(name, msa_rows, spans=spans)

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
            pae_rows = format_pae_rows(data["pae"])
            token_rows = format_token_chain_rows(data.get("token_chain_ids"), data["pae"])
            write_tsv(f"{name}_{rank_id}_pae.tsv", token_rows + pae_rows)
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
    parser.add_argument(
        "--af3_msa_json", dest="af3_msa_json", required=False,
        help="AlphaFold3-format merged MSA JSON; expands entity copies and "
             "orders the MSA columns by chain id against --structure",
    )
    parser.add_argument(
        "--structure", dest="structure", required=False,
        help="Top-ranked predicted structure; orders the expanded MSA columns "
             "into the structure's polymer-chain order so the report's MSA "
             "coverage aligns with the PAE tokens",
    )
    parser.add_argument("--colabfold_metrics_files", required=False, nargs="+")
    parser.add_argument("--structs", dest="structs", required=False, nargs="+")
    parser.add_argument("--name", default="untitled", dest="name") # might need a --name $meta.id
    args = parser.parse_args()

    if args.af3_msa_json:
        read_af3_msa_json(args.name, args.af3_msa_json, structure=args.structure)
    elif args.paired_a3m:
        read_paired_a3m(args.name, args.paired_a3m, structure=args.structure)
    if args.structs:
        extract_structs_plddt_to_tsv(args.name, args.structs, args.colabfold_metrics_files)
    if args.colabfold_metrics_files:
        read_colabfold_metrics(args.name, args.colabfold_metrics_files)

if __name__ == "__main__":
    main()
