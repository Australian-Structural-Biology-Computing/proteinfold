#!/usr/bin/env python3
"""Convert AlphaFold2 MSA pipeline output into an AlphaFold3-format input JSON.

The unified proteinfold workflow runs structure prediction exclusively through
ColabFold2 (colabfold_batch). colabfold_batch accepts an AlphaFold3-format job
JSON whose protein entities carry `unpairedMsa` / `pairedMsa` a3m text blocks
(this is what `colabfold_search --af3-json` produces). This script rebuilds
that JSON from the output of RUN_ALPHAFOLD2_MSA (alphafold's run_msa.py):

  <msa_dir>/chain_id_map.json
  <msa_dir>/<chain>/uniref90_hits.sto      (jackhmmer, unclustered uniref90)
  <msa_dir>/<chain>/bfd_uniref_hits.a3m    (hhblits bfd+uniref30, full dbs)
    or <msa_dir>/<chain>/small_bfd_hits.sto (reduced dbs)
  <msa_dir>/<chain>/mgnify_hits.sto        (jackhmmer mgnify)
  <msa_dir>/<chain>/uniprot_hits.sto       (jackhmmer uniprot, pairing only,
                                            present when the input has >= 2
                                            unique sequences)

unpairedMsa: the chain's own MSA, assembled exactly like alphafold's
  pipeline.DataPipeline: uniref90 -> bfd -> mgnify, deduplicated by aligned
  sequence, emitted as a3m text (query row first, lowercase insertions).

pairedMsa: alphafold's pipeline_multimer species pairing (msa_pairing.py
  pair_sequences / reorder_paired_rows), ported here, applied to the uniprot
  jackhmmer hits. Output rows are re-described with synthetic UniProt-style
  accessions so that row *i* of every chain pairs together (the convention
  colabfold/msa_pairing.py establishes for alphafold3 species pairing).
  For homomers (all chains identical) alphafold skips the uniprot search, so
  the uniref90-based unpaired MSA is reused as the paired block: row *i* of
  identical chains necessarily corresponds to the same species.
  For monomers pairedMsa is left empty.

Stockholm inputs compressed with zstd (*.sto.zst) are transparently
decompressed via the `zstd` binary.
"""

import argparse
import json
import os
import re
import subprocess
import sys

# Same caps the MMseqs search applies to colabfold_search output.
MAX_MSA_SEQS = 16384
MAX_PAIRED_SEQS = 8192

# alphafold/data/msa_pairing.py: skip species with more than this many hits.
MAX_PER_SPECIES_ROWS = 600

_UNIPROT_PATTERN = re.compile(
    r"""
    ^
    (?:tr|sp)
    \|(?P<AccessionIdentifier>[A-Za-z0-9]{6,10})(?:_\d)?
    \|(?:[A-Za-z0-9]+)_(?P<SpeciesIdentifier>[A-Za-z0-9]{1,5})(?:_\d+)?
    $
    """,
    re.VERBOSE,
)


def parse_fasta(fasta_string):
    """Port of alphafold.data.parsers.parse_fasta."""
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        if not line:
            continue
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])
            sequences.append("")
            continue
        if index < 0:
            raise ValueError(f"Sequence data before FASTA header: {line}")
        sequences[index] += line
    return sequences, descriptions


def read_text(path):
    """Read a possibly zstd-compressed text file."""
    if path.endswith(".zst"):
        result = subprocess.run(
            ["zstd", "-dc", path], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        return result.stdout.decode("utf-8")
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def find_alignment_file(directory, name):
    """Locate an alignment file, accepting a .zst-compressed variant.

    Institutional RUN_ALPHAFOLD2_MSA configs may compress the Stockholm
    outputs with zstd.
    """
    for candidate in (name, name + ".zst"):
        path = os.path.join(directory, candidate)
        if os.path.exists(path):
            return path
    return None


def stockholm_to_rows(stockholm_string):
    """Port of parsers.convert_stockholm_to_a3m, returning (desc, seq) rows.

    The first row is the query; gaps present in the query row are removed and
    residues aligned there in other rows become lowercase insertions.
    """
    sequences = {}
    descriptions = {}
    for line in stockholm_string.splitlines():
        if line.strip() and not line.startswith(("#", "//")):
            seqname, aligned_seq = line.split(maxsplit=1)
            if seqname not in sequences:
                sequences[seqname] = ""
            sequences[seqname] += aligned_seq
    for line in stockholm_string.splitlines():
        if line[:4] == "#=GS":
            columns = line.split(maxsplit=3)
            seqname, feature = columns[1:3]
            if feature != "DE":
                continue
            descriptions[seqname] = columns[3] if len(columns) == 4 else ""

    if not sequences:
        return []

    query_sequence = next(iter(sequences.values()))
    query_non_gaps = [res != "-" for res in query_sequence]

    rows = []
    for seqname, sto_sequence in sequences.items():
        out_sequence = sto_sequence.replace(".", "")
        converted = []
        for is_query_res_non_gap, sequence_res in zip(query_non_gaps, out_sequence):
            if is_query_res_non_gap:
                converted.append(sequence_res)
            elif sequence_res != "-":
                converted.append(sequence_res.lower())
        description = descriptions.get(seqname, "")
        rows.append((f"{seqname} {description}".strip(), "".join(converted)))
    return rows


def a3m_to_rows(a3m_string):
    """Parse an a3m string into (description, sequence) rows."""
    rows = []
    for line in a3m_string.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            rows.append([line[1:], ""])
        elif rows:
            rows[-1][1] += line
    return [tuple(row) for row in rows]


def aligned_sequence(row_sequence):
    """Deletion-free, gap-carrying representation used for dedup/pairing."""
    return "".join(
        residue
        for residue in row_sequence
        if not residue.islower()
    )


def species_identifier(description):
    """Port of alphafold.data.msa_identifiers.get_identifiers."""
    split_description = description.split()
    if not split_description:
        return ""
    sequence_identifier = split_description[0].partition("/")[0]
    matches = _UNIPROT_PATTERN.match(sequence_identifier.strip())
    if matches:
        return matches.group("SpeciesIdentifier")
    return ""


def rows_to_a3m(rows):
    return "\n".join(
        f">{description}\n{sequence}" for description, sequence in rows
    ) + "\n"


def deduplicate_rows(rows):
    """Port of pipeline.make_msa_features dedup: first occurrence wins."""
    seen = set()
    output = []
    for description, sequence in rows:
        key = aligned_sequence(sequence)
        if key in seen:
            continue
        seen.add(key)
        output.append((description, sequence))
    return output


def build_unpaired_rows(chain_msa_dir, use_small_bfd):
    """Assemble the chain's unpaired MSA: uniref90 -> bfd -> mgnify."""
    paths = {
        "uniref90": find_alignment_file(chain_msa_dir, "uniref90_hits.sto"),
        "bfd_a3m": find_alignment_file(chain_msa_dir, "bfd_uniref_hits.a3m"),
        "small_bfd": find_alignment_file(chain_msa_dir, "small_bfd_hits.sto"),
        "mgnify": find_alignment_file(chain_msa_dir, "mgnify_hits.sto"),
    }
    rows = []
    for name in ("uniref90", "bfd_a3m", "small_bfd", "mgnify"):
        path = paths[name]
        if name == "bfd_a3m":
            if use_small_bfd or path is None:
                continue
            rows.extend(a3m_to_rows(read_text(path)))
        else:
            if (name == "small_bfd") != use_small_bfd:
                continue
            if path is None:
                if name in ("uniref90", "mgnify"):
                    raise FileNotFoundError(
                        f"Missing expected alignment file: "
                        f"{os.path.join(chain_msa_dir, name)}"
                    )
                continue
            rows.extend(stockholm_to_rows(read_text(path)))
    if not rows:
        raise ValueError(f"No alignment rows found in {chain_msa_dir}")
    return deduplicate_rows(rows)


def build_paired_rows(chain_msa_dir):
    """Parse the chain's uniprot jackhmmer hits for species pairing."""
    path = find_alignment_file(chain_msa_dir, "uniprot_hits.sto")
    if path is None:
        return None
    return stockholm_to_rows(read_text(path))


def pair_species_chains(chain_paired_rows):
    """Port of alphafold.data.msa_pairing.pair_sequences.

    Args:
        chain_paired_rows: list (one entry per chain) of lists of
            (description, a3m_sequence) rows from uniprot_hits.sto.

    Returns:
        A list of paired rows; each row is a list of per-chain indices into
        that chain's rows, where an index equal to the chain's row count
        denotes the padding (all-gap) row. Row 0 of the result is the query.
    """
    num_chains = len(chain_paired_rows)

    chain_species = []
    common_species = set()
    for rows in chain_paired_rows:
        query = aligned_sequence(rows[0][1])
        query_length = len(query)
        species_rows = {}
        for row_index, (description, sequence) in enumerate(rows):
            aligned = aligned_sequence(sequence)
            similarity = sum(
                query_res == row_res
                for query_res, row_res in zip(query, aligned)
            ) / query_length
            species = species_identifier(description)
            species_rows.setdefault(species, []).append(
                (row_index, similarity)
            )
        chain_species.append(species_rows)
        common_species.update(species_rows)

    common_species = sorted(common_species)
    if "" in common_species:
        common_species.remove("")  # The query sequence itself.

    paired_rows = [[0] * num_chains]
    rows_by_chain_count = {num_chains: [tuple([0] * num_chains)]}

    for species in common_species:
        if not species:
            continue
        present = 0
        species_dfs = []
        for species_rows in chain_species:
            if species in species_rows:
                species_dfs.append(species_rows[species])
                present += 1
            else:
                species_dfs.append(None)

        # Skip species present in only one chain.
        if present <= 1:
            continue
        # Skip species with an excessive number of hits in any chain.
        if any(len(rows) > MAX_PER_SPECIES_ROWS for rows in species_dfs if rows):
            continue

        take_count = min(len(rows) for rows in species_dfs if rows)
        row_indices = []
        for rows in species_dfs:
            if rows is None:
                # Take the padding row.
                row_indices.append([None] * take_count)
            else:
                # Most-similar-to-query first; ties keep input order.
                sorted_rows = sorted(
                    rows, key=lambda item: -item[1]
                )[:take_count]
                row_indices.append([item[0] for item in sorted_rows])
        for indices in zip(*row_indices):
            paired_rows.append(list(indices))
            rows_by_chain_count.setdefault(present, []).append(tuple(indices))

    # Port of reorder_paired_rows: group rows by how many chains they pair,
    # most chains first; within a group order by product of row indices.
    ordered = []
    for chain_count in sorted(rows_by_chain_count, reverse=True):
        group = rows_by_chain_count[chain_count]
        group.sort(key=lambda indices: _row_product(indices, chain_paired_rows))
        ordered.extend(group)
    return ordered


def _row_product(indices, chain_paired_rows):
    """abs(product of row indices), padding rows use their appended index."""
    product = 1
    for chain_num, index in enumerate(indices):
        if index is None:
            index = len(chain_paired_rows[chain_num])
        product *= index
    return abs(product)


def paired_blocks_from_indices(paired_indices, chain_paired_rows):
    """Materialise per-chain paired a3m blocks from pairing row indices."""
    blocks = []
    for chain_num, rows in enumerate(chain_paired_rows):
        lines = [f">{101}"]  # Row 0 is the query; empty species id.
        lines.append(rows[0][1])
        block_row = 0
        for indices in paired_indices[1:]:
            index = indices[chain_num]
            if index is None:
                sequence = "-" * len(aligned_sequence(rows[0][1]))
            else:
                sequence = rows[index][1]
            block_row += 1
            # Synthetic description: row *i* pairs across all chains.
            lines.append(f">tr|AAAAAA|AAAAAA_{block_row:05d}")
            lines.append(sequence)
        blocks.append("\n".join(lines) + "\n")
    return blocks


def build_job(name, chains, chain_ids, unpaired, paired):
    """Assemble the AlphaFold3-format job JSON (colabfold AF3Utils layout)."""
    sequences = []
    for chain_id, sequence, unpaired_msa, paired_msa in zip(
        chain_ids, chains, unpaired, paired
    ):
        protein = {
            "id": [chain_id],
            "sequence": sequence,
            "modifications": [],
            "templates": [],
        }
        protein["unpairedMsa"] = unpaired_msa if unpaired_msa else ""
        protein["pairedMsa"] = paired_msa if paired_msa else ""
        sequences.append({"protein": protein})
    return {
        "dialect": "alphafold3",
        "version": 2,
        "name": name,
        "sequences": sequences,
        "modelSeeds": [1],
        "bondedAtomPairs": None,
        "userCCD": None,
    }


def int_id_to_str_id(number):
    """Port of alphafold.int_id_to_str_id (A, B, ..., Z, AA, ...)."""
    if number <= 0:
        raise ValueError(f"Only positive integers allowed, got {number}.")
    number -= 1
    output = []
    while number >= 0:
        output.append(chr(number % 26 + ord("A")))
        number = number // 26 - 1
    return "".join(output)


def discover_msa_dir(path):
    """Locate the MSA output directory under a staged dir.

    RUN_ALPHAFOLD2_MSA writes <output>/<fasta_basename>/msas/, but staging may
    flatten or rename levels, so accept the dir itself, msas/, or any unique
    match up to three levels deep. The multimer data pipeline writes a
    chain_id_map.json; the monomer pipeline writes the alignments at the MSA
    directory root instead. Returns (msa_dir, is_multimer_layout).
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
        return map_candidates[0], True
    if not map_candidates and len(root_candidates) == 1:
        return root_candidates[0], False
    raise ValueError(
        f"Expected exactly one MSA output directory under {path} "
        f"(chain_id_map.json or uniref90_hits.sto), found "
        f"chain_id_map.json: {map_candidates}, uniref90: {root_candidates}"
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument(
        "--msa_dir",
        required=True,
        help="Directory holding the run_msa.py output for this FASTA "
        "(contains chain_id_map.json, either directly or under msas/).",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--name", required=True, help="Job/sample name.")
    parser.add_argument(
        "--max-unpaired", type=int, default=MAX_MSA_SEQS
    )
    parser.add_argument(
        "--max-paired", type=int, default=MAX_PAIRED_SEQS
    )
    parser.add_argument(
        "--small-bfd", action="store_true",
        help="RUN_ALPHAFOLD2_MSA ran with the reduced_dbs preset.",
    )
    args = parser.parse_args()

    sequences, _ = parse_fasta(open(args.fasta, encoding="utf-8").read())

    msa_dir, multimer_layout = discover_msa_dir(args.msa_dir)

    if multimer_layout:
        with open(os.path.join(msa_dir, "chain_id_map.json"), encoding="utf-8") as handle:
            chain_id_map = json.load(handle)
        chain_ids = sorted(chain_id_map)
        if len(chain_ids) != len(sequences):
            raise ValueError(
                f"FASTA has {len(sequences)} chain(s) but chain_id_map.json has "
                f"{len(chain_ids)}: {chain_ids}"
            )
        # Chains with an identical sequence share one MSA directory.
        chain_dirs = {}
        for chain_id, fasta_chain in chain_id_map.items():
            chain_dirs.setdefault(fasta_chain["sequence"], chain_id)
        chain_sequences = {
            chain_id: chain_id_map[chain_id]["sequence"] for chain_id in chain_ids
        }
    else:
        # Monomer data pipeline: single chain, alignments at the MSA dir root.
        if len(sequences) != 1:
            raise ValueError(
                "The monomer data pipeline ran for a multi-chain entry; "
                f"{len(sequences)} chains cannot be mapped"
            )
        chain_ids = ["A"]
        chain_dirs = {sequences[0]: "A"}
        chain_sequences = {"A": sequences[0]}

    is_homomer_or_monomer = len(set(sequences)) == 1
    unique_chain_dirs = {
        chain_id_map[chain_id]["sequence"]: chain_id_map[chain_id]
        for chain_id in chain_ids
    } if multimer_layout else {}
    has_pairing_input = (
        not is_homomer_or_monomer and len(unique_chain_dirs) > 1
    )

    unpaired_rows_by_chain = {}
    paired_rows_by_chain = {}
    for chain_id in chain_ids:
        sequence = chain_sequences[chain_id]
        shared_chain = chain_dirs[sequence]
        chain_msa_dir = (
            os.path.join(msa_dir, shared_chain) if multimer_layout else msa_dir
        )
        if chain_id not in unpaired_rows_by_chain:
            unpaired_rows_by_chain[chain_id] = build_unpaired_rows(
                chain_msa_dir, args.small_bfd
            )
            paired_rows_by_chain[chain_id] = build_paired_rows(chain_msa_dir)
        else:
            unpaired_rows_by_chain[chain_id] = unpaired_rows_by_chain[
                shared_chain
            ]
            paired_rows_by_chain[chain_id] = paired_rows_by_chain[shared_chain]

    # Cap the unpaired blocks, always keeping the query row.
    unpaired_blocks = []
    for chain_id in chain_ids:
        rows = unpaired_rows_by_chain[chain_id]
        kept = rows[:1] + rows[1 : args.max_unpaired]
        unpaired_blocks.append(rows_to_a3m(kept))

    paired_blocks = None
    if has_pairing_input:
        ordered_paired = []
        for chain_id in chain_ids:
            rows = paired_rows_by_chain[chain_id]
            if rows is None:
                raise ValueError(
                    f"Chain {chain_id} has no uniprot_hits.sto; MSA pairing "
                    "requires a multimer-preset RUN_ALPHAFOLD2_MSA run."
                )
            ordered_paired.append(rows)
        paired_indices = pair_species_chains(ordered_paired)
        if len(paired_indices) - 1 > args.max_paired:
            print(
                f"Truncating paired MSA from {len(paired_indices) - 1} to "
                f"{args.max_paired} rows",
                file=sys.stderr,
            )
            paired_indices = paired_indices[: args.max_paired + 1]
        paired_blocks = paired_blocks_from_indices(
            paired_indices, ordered_paired
        )
    elif is_homomer_or_monomer and len(sequences) > 1:
        # Homomer: identical chains, so row *i* of each chain's unpaired MSA
        # corresponds to the same species. Reuse it as the paired block.
        depth = min(
            len(unpaired_rows_by_chain[chain_id]) for chain_id in chain_ids
        )
        depth = min(depth, args.max_paired)
        paired_blocks = [
            rows_to_a3m(unpaired_rows_by_chain[chain_id][:depth])
            for chain_id in chain_ids
        ]

    job = build_job(
        args.name,
        [chain_sequences[chain_id] for chain_id in chain_ids],
        [int_id_to_str_id(index + 1) for index in range(len(chain_ids))],
        unpaired_blocks,
        paired_blocks or [""] * len(chain_ids),
    )

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(job, indent=2) + "\n")


if __name__ == "__main__":
    main()
