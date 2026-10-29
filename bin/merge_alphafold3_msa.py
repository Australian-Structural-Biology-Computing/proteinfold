#!/usr/bin/env python3

import argparse
import copy
import json
from pathlib import Path


def _single_job(document, label):
    if isinstance(document, list):
        if len(document) != 1:
            raise ValueError(f"{label} must contain exactly one AlphaFold3 job")
        document = document[0]
    if not isinstance(document, dict):
        raise ValueError(f"{label} must be an AlphaFold3 JSON object")
    return document


def _proteins(document):
    proteins = []
    for entity in document.get("sequences", []):
        if isinstance(entity, dict) and isinstance(entity.get("protein"), dict):
            proteins.append(entity["protein"])
    return proteins


def merge_msa(original, searched):
    merged = copy.deepcopy(original)
    available = _proteins(searched)

    for protein in _proteins(merged):
        sequence = protein.get("sequence")
        candidate = next(
            (
                item
                for item in available
                if str(item.get("sequence", "")).upper() == str(sequence or "").upper()
            ),
            None,
        )
        if candidate is None:
            raise ValueError(
                "MMseqs output has no protein entity matching sequence "
                f"for original chain(s) {protein.get('id')}"
            )
        for field in ("unpairedMsa", "pairedMsa"):
            if protein.get(field) is not None:
                continue
            value = candidate.get(field)
            if value is None:
                raise ValueError(
                    f"MMseqs protein entity {candidate.get('id')} is missing {field}"
                )
            protein[field] = value

    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Merge MMseqs-generated MSAs into an AlphaFold3 input JSON"
    )
    parser.add_argument("msa_json")
    parser.add_argument("original_json")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.msa_json, encoding="utf-8") as handle:
        searched = _single_job(json.load(handle), "MMseqs JSON")
    with open(args.original_json, encoding="utf-8") as handle:
        original = _single_job(json.load(handle), "original JSON")

    output = merge_msa(original, searched)
    Path(args.output).write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
