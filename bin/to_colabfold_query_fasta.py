#!/usr/bin/env python3
"""Reduce a sample's input to the single-query ColabFold FASTA.

colabfold_search treats a sequence containing ':' as a complex, so every
chain of a sample is colon-joined into one query record: multi-record
FASTAs are joined across their records, and AlphaFold3 JSONs are reduced
to their protein chains (repeated per copy, in entity order), so the
search sees the same chains the JSON describes. Non-protein entities
(RNA, DNA, ligands) are ignored here; for JSON inputs the original
document is repopulated with the generated MSAs downstream
(merge_alphafold3_msa.py).
"""

import argparse
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


def fasta_sequences(path):
    """Read a multi-record FASTA as one sequence per record."""
    sequences = []
    current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current is not None:
                sequences.append("".join(current))
            current = []
        elif line:
            if current is None:
                raise ValueError(f"FASTA input {path} has sequence data before any header line")
            current.append(line)
    if current is not None:
        sequences.append("".join(current))
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description="Reduce a sample FASTA or AlphaFold3 JSON to the single-query ColabFold FASTA"
    )
    parser.add_argument("input")
    parser.add_argument("--id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if input_path.suffix == ".json":
        queries = protein_queries(_single_job(json.loads(input_path.read_text(encoding="utf-8"))))
    else:
        queries = fasta_sequences(input_path)
        if not queries:
            raise ValueError(f"FASTA input {input_path} contains no records")
    sequence = ":".join(queries)
    Path(args.output).write_text(f">{args.id}\n{sequence}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
