#!/usr/bin/env python3
"""Extract the protein chains of an AlphaFold3 JSON into a multi-record FASTA.

The unified proteinfold pipeline's alphafold2_msa MSA tool runs alphafold's
data pipeline, which expects one FASTA record per chain (a multi-record FASTA
is treated as a multimer and species-paired across its chains). This script
expands every protein entity of a single AlphaFold3 job into one record per
copy, in entity order, so the MSA pipeline sees the same chains the JSON
describes. Non-protein entities (RNA, DNA, ligands) are ignored here; the
original JSON is repopulated with the generated MSAs downstream
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


def protein_chains(document):
    chains = []
    for entity in document.get("sequences", []):
        if not isinstance(entity, dict) or "protein" not in entity:
            continue
        protein = entity["protein"]
        sequence = protein.get("sequence") if isinstance(protein, dict) else None
        if not sequence:
            raise ValueError("AlphaFold3 protein entity is missing its sequence")
        entity_ids = protein.get("id")
        copies = len(entity_ids) if isinstance(entity_ids, list) else 1
        chains.extend([sequence] * max(1, copies))
    if not chains:
        raise ValueError(
            "alphafold2_msa MSA generation needs at least one protein entity; "
            "use --msa_tool colabfold_search or alphafold3_msa for an "
            "RNA/ligand-only AlphaFold3 JSON"
        )
    return chains


def main():
    parser = argparse.ArgumentParser(
        description="Convert the protein entities in AlphaFold3 JSON to a multi-record FASTA for the AlphaFold2 MSA pipeline"
    )
    parser.add_argument("json_in")
    parser.add_argument("--id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.json_in, encoding="utf-8") as handle:
        document = _single_job(json.load(handle))

    chunks = []
    for index, sequence in enumerate(protein_chains(document), start=1):
        chunks.append(f">{args.id}_{index}\n{sequence}")
    Path(args.output).write_text("\n".join(chunks) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
