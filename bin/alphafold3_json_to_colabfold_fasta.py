#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


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


def protein_queries(document):
    queries = []
    for entity in document.get("sequences", []):
        if not isinstance(entity, dict) or "protein" not in entity:
            continue
        protein = entity["protein"]
        sequence = protein.get("sequence") if isinstance(protein, dict) else None
        if not sequence:
            raise ValueError("AlphaFold3 protein entity is missing its sequence")
        entity_ids = protein.get("id")
        copies = len(entity_ids) if isinstance(entity_ids, list) else 1
        queries.extend([sequence] * max(1, copies))
    if not queries:
        raise ValueError(
            "Offline ColabFold2 MSA generation needs at least one protein entity; "
            "use --use_msa_server for an RNA/ligand-only AlphaFold3 JSON"
        )
    return queries


def main():
    parser = argparse.ArgumentParser(
        description="Convert the protein entities in AlphaFold3 JSON to a ColabFold query CSV"
    )
    parser.add_argument("json_in")
    parser.add_argument("--id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.json_in, encoding="utf-8") as handle:
        document = _single_job(json.load(handle))
    sequence = ":".join(protein_queries(document))
    # The MMseqs search (and the MSA cache) consume id,sequence query CSVs;
    # colon-joined chains keep complexes searched (and paired) as one query.
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "sequence"])
        writer.writerow([args.id, sequence])


if __name__ == "__main__":
    main()
