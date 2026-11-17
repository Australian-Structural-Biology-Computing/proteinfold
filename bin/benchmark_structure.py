#!/usr/bin/env python3
"""Compare a single-chain prediction with an experimental structure."""

import argparse
import json
import math
import subprocess
import tempfile
from pathlib import Path

from ost import io


def check_protein(path):
    structure = io.LoadPDB(str(path)) if path.suffix.lower() == ".pdb" else io.LoadMMCIF(str(path))
    peptides = structure.Select("peptide=true")
    nucleotides = structure.Select("nucleotide=true")
    if len(peptides.chains) != 1 or len(nucleotides.residues):
        raise ValueError(f"{path}: benchmarking requires exactly one protein chain and no nucleic acid chains")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    check_protein(args.model)
    check_protein(args.reference)
    with tempfile.TemporaryDirectory() as tmp:
        result = Path(tmp) / "scores.json"
        subprocess.run(
            [
                "ost", "compare-structures", "-m", str(args.model), "-r", str(args.reference),
                "--lddt", "--tm-score", "--enforce-consistency", "-o", str(result),
            ],
            check=True,
        )
        data = json.loads(result.read_text())
    if data.get("status") != "SUCCESS":
        raise ValueError(f"OpenStructure comparison failed: {data.get('traceback', data)}")
    if len(data["reference_chains"]) != 1 or len(data["model_chains"]) != 1 or len(data["chain_mapping"]) != 1:
        raise ValueError("OpenStructure could not uniquely map the single protein chains")
    if data.get("usalign_mapping") != data["chain_mapping"]:
        raise ValueError("OpenStructure and US-align produced different chain mappings")
    for metric in ("lddt", "tm_score"):
        score = data.get(metric)
        if not isinstance(score, (int, float)) or not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError(f"OpenStructure returned an invalid {metric}: {score!r}")
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
