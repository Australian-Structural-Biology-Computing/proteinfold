#!/usr/bin/env python3
"""Convert an AlphaFold3-format MSA JSON into a ColabFold a3m query file.

The colabfold2 alphafold2 backends (colabfold_batch --model-type
alphafold2*) predate the AlphaFold3 input JSON: they consume a ColabFold
a3m file (the format the ColabFold webserver and colabfold_search emit) and
ignore any MSAs embedded in an AlphaFold3 JSON. The unified proteinfold
pipeline generates MSAs with one of three tools, all of which publish an
AlphaFold3-format JSON; this adaptor re-serialises those MSAs into the
combined a3m format so the alphafold2 backends can consume them.

The output format replicates colabfold.input.msa_to_str:

    #<len1>,...,<lenN>\t<card1>,...,<cardN>
    >101\t101
    <seq1><seq2>...          (paired block; tab-joined chain descriptions)
    ...
    >101
    <chain 1 unpaired row>   (unpaired block; '-'-padded to full width)
    ...

Non-protein entities (RNA, ligands, ...) are skipped: the alphafold2
backends fold protein chains only.
"""

import argparse
import json
import sys


def _single_job(document):
    if isinstance(document, list):
        if len(document) != 1:
            raise ValueError(
                "ColabFold2 requires one AlphaFold3 job per samplesheet row; "
                f"the JSON contains {len(document)} jobs"
            )
        document = document[0]
    if not isinstance(document, dict):
        raise ValueError("AlphaFold3 JSON must be an object or a one-element list")
    return document


def protein_chains(document):
    """Unique protein sequences (first-appearance order) and their copy counts."""
    unique = []
    counts = []
    for entity in document.get("sequences", []):
        if not isinstance(entity, dict) or "protein" not in entity:
            continue
        protein = entity["protein"]
        sequence = protein.get("sequence") if isinstance(protein, dict) else None
        if not sequence:
            raise ValueError("AlphaFold3 protein entity is missing its sequence")
        entity_ids = protein.get("id")
        copies = len(entity_ids) if isinstance(entity_ids, list) else 1
        copies = max(1, copies)
        if sequence in unique:
            counts[unique.index(sequence)] += copies
        else:
            unique.append(sequence)
            counts.append(copies)
    if not unique:
        raise ValueError(
            "the AlphaFold2 backends fold protein chains only; the JSON has no "
            "protein entity"
        )
    return unique, counts


def a3m_rows(a3m_text):
    """Parse an a3m text block into (description, sequence) rows."""
    rows = []
    for line in (a3m_text or "").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            rows.append([line[1:], ""])
        elif rows:
            rows[-1][1] += line
    return rows


def pair_sequences(paired_msas, cardinalities):
    """Port of colabfold.input.pair_sequences (row-order pairing)."""
    depths = [len(a3m.splitlines()) for a3m in paired_msas]
    if len(set(depths)) != 1:
        raise ValueError(
            "cannot pair MSA blocks of unequal depth: {}; the AlphaFold2 "
            "backends require row-aligned paired MSAs".format(depths)
        )
    lines = [""] * depths[0]
    for chain, a3m in enumerate(paired_msas):
        for index, line in enumerate(a3m.splitlines()):
            if line.startswith(">"):
                if chain != 0:
                    line = line.replace(">", "\t", 1)
                lines[index] = lines[index] + line
            else:
                lines[index] = lines[index] + line * cardinalities[chain]
    return "\n".join(lines)


def pad_sequences(unpaired_msas, sequences, cardinalities):
    """Port of colabfold.input.pad_sequences ('-'-padding across chains)."""
    blank_seq = []
    for position, seq in enumerate(sequences):
        for _ in range(cardinalities[position]):
            blank_seq.append("-" * len(seq))

    combined = []
    position = 0
    for chain, seq in enumerate(sequences):
        for _ in range(cardinalities[chain]):
            for line in unpaired_msas[chain].split("\n"):
                if len(line) == 0:
                    continue
                if line.startswith(">"):
                    combined.append(line)
                else:
                    combined.append(
                        "".join(blank_seq[:position] + [line] + blank_seq[position + 1:])
                    )
            position += 1
    return "\n".join(combined)


def main():
    parser = argparse.ArgumentParser(
        description="Convert an AlphaFold3-format MSA JSON into a ColabFold a3m query file"
    )
    parser.add_argument("json_in")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.json_in, encoding="utf-8") as handle:
        document = _single_job(json.load(handle))

    sequences, cardinalities = protein_chains(document)

    unpaired_msas = []
    paired_msas = []
    have_paired = False
    for entity in document.get("sequences", []):
        if not isinstance(entity, dict) or "protein" not in entity:
            continue
        protein = entity["protein"]
        sequence = protein.get("sequence")
        if sequence not in sequences:
            continue  # only the first copy of each unique sequence carries the MSA
        unpaired = protein.get("unpairedMsa")
        paired = protein.get("pairedMsa")
        if not unpaired or not unpaired.strip():
            raise ValueError(
                "protein chain {} has no unpairedMsa; the AlphaFold2 backends "
                "need a full unpaired MSA".format(protein.get("id"))
            )
        unpaired_msas.append(unpaired)
        if paired and paired.strip():
            paired_msas.append(paired)
            have_paired = True
        else:
            paired_msas.append("")

    header = "#{}\t{}".format(
        ",".join(str(len(seq)) for seq in sequences),
        ",".join(str(count) for count in cardinalities),
    )

    if have_paired:
        blank_chains = [i for i, a3m in enumerate(paired_msas) if not a3m.strip()]
        if blank_chains:
            raise ValueError(
                "protein chains {} have no pairedMsa while others do; the "
                "AlphaFold2 backends cannot mix paired and unpaired-only "
                "chains - re-run the MSA search with pairing enabled".format(
                    [sequences[i] for i in blank_chains]
                )
            )
        body = pair_sequences(paired_msas, [1] * len(sequences)) + "\n" + pad_sequences(
            unpaired_msas, sequences, [1] * len(sequences)
        )
    else:
        body = pad_sequences(unpaired_msas, sequences, [1] * len(sequences))

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(header + "\n" + body + "\n")


if __name__ == "__main__":
    main()
