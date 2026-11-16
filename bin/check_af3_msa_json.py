#!/usr/bin/env python3
"""Validate an AlphaFold3-format JSON produced by a pre-computed MSA pipeline.

The unified proteinfold workflow feeds this JSON straight to ColabFold2
(colabfold_batch). colabfold reads `unpairedMsa` / `pairedMsa` from each
protein chain: a string of a3m content, "" meaning deliberately empty and a
missing/None field meaning "run the MSA search". Since the point of the
alphafold3_msa MSA tool is to supply pre-computed MSAs, this script refuses
JSONs where a protein chain would silently trigger an extra search.
"""

import argparse
import json
import sys


def _single_job(document, label):
    if isinstance(document, list):
        if len(document) != 1:
            raise ValueError(f"{label} must contain exactly one AlphaFold3 job")
        document = document[0]
    if not isinstance(document, dict):
        raise ValueError(f"{label} must be an AlphaFold3 JSON object")
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.input_json, encoding="utf-8") as handle:
        document = _single_job(json.load(handle), "input JSON")

    problems = []
    proteins = 0
    for entity in document.get("sequences", []):
        protein = entity.get("protein") if isinstance(entity, dict) else None
        if not isinstance(protein, dict):
            continue
        proteins += 1
        chain_ids = protein.get("id")
        label = f"protein chain {chain_ids}"
        unpaired = protein.get("unpairedMsa")
        paired = protein.get("pairedMsa")
        if unpaired is None and paired is None:
            problems.append(
                f"{label} has neither unpairedMsa nor pairedMsa; "
                "ColabFold2 would run its own MSA search"
            )
        elif unpaired is None:
            problems.append(
                f"{label} has a pairedMsa but no unpairedMsa; ColabFold2 "
                "would run an extra unpaired MSA search"
            )
        for field, value in (("unpairedMsa", unpaired), ("pairedMsa", paired)):
            if value is not None and not isinstance(value, str):
                problems.append(f"{label} {field} must be an a3m text string")
        sequence = protein.get("sequence")
        if not sequence:
            problems.append(f"{label} has no sequence")

    if proteins == 0:
        problems.append("JSON contains no protein chains")

    if problems:
        for problem in problems:
            print(f"ERROR: {problem}", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(document, indent=2) + "\n")


if __name__ == "__main__":
    main()
