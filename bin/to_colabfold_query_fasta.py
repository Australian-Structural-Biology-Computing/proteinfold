#!/usr/bin/env python3
"""Reduce sample inputs to colon-joined ColabFold query FASTAs.

colabfold_search treats a sequence containing ':' as a complex, so every
chain of a sample is colon-joined into one query record: multi-record
FASTAs are joined across their records, and AlphaFold3 JSONs are reduced
to their protein chains (repeated per copy, in entity order), so the
search sees the same chains the JSON describes. Non-protein entities
(RNA, DNA, ligands) are ignored here; for JSON inputs the original
document is repopulated with the generated MSAs downstream
(merge_alphafold3_msa.py).

Subcommands:
  single INPUT --id ID --output PATH
      Reduce one sample input to its single-query FASTA.
  batch --manifest TSV [--batch-size N] (--chunks-dir D | --queries-dir D)
      Process a whole run in one go without per-sample intermediates:
      the queries are concatenated on the fly into chunk FASTAs of
      --batch-size records, or written one per sample with --queries-dir.
"""


import argparse
import json
import os
import sys
from pathlib import Path


# --------------------------------------------------------------------------
# Query conversion
# --------------------------------------------------------------------------


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


def sample_queries(input_path):
    """Colon-joined query sequence for one sample input (FASTA or AF3 JSON)."""
    input_path = Path(input_path)
    if input_path.suffix == ".json":
        queries = protein_queries(_single_job(json.loads(input_path.read_text(encoding="utf-8"))))
    else:
        queries = fasta_sequences(input_path)
        if not queries:
            raise ValueError(f"FASTA input {input_path} contains no records")
    return ":".join(queries)


def _staged_path(name):
    """Task-context paths may already carry the input/ staging prefix."""
    return name if name.startswith("input/") else "input/" + name


def _read_manifest(path):
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(f"Manifest line must be '<sample id>\\t<record header>\\t<input file>': {line!r}")
        rows.append(tuple(fields))
    if not rows:
        raise ValueError(f"Manifest {path} contains no samples")
    return rows


def run_batch(args):
    if args.batch_size == 0 or args.batch_size < -1:
        raise ValueError(f"--batch-size must be -1 or a positive integer, got {args.batch_size}")
    if bool(args.queries_dir) == bool(args.chunks_dir):
        raise ValueError("Choose exactly one of --queries-dir and --chunks-dir")
    rows = _read_manifest(args.manifest)
    sample_ids = [row[0] for row in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError(f"Manifest contains duplicate sample ids: {sample_ids}")
    for row in rows:
        for field in row:
            if any(character in field for character in "\t\r\n"):
                raise ValueError(f"Manifest fields must not contain tabs or newlines: {row}")

    # Deterministic chunk composition regardless of channel arrival order.
    rows = sorted(rows, key=lambda row: row[0])

    entries = []
    for sample_id, record_id, input_name in rows:
        entries.append((sample_id, record_id, sample_queries(_staged_path(input_name))))

    queries = entries

    if args.queries_dir:
        os.makedirs(args.queries_dir, exist_ok=True)
        for sample_id, record_id, sequence in queries:
            Path(os.path.join(args.queries_dir, sample_id + ".fasta")).write_text(
                f">{record_id}\n{sequence}\n", encoding="utf-8"
            )
        return

    limit = args.batch_size if args.batch_size > 0 else None
    index = 0
    buffer = []

    def flush():
        os.makedirs(args.chunks_dir, exist_ok=True)
        Path(os.path.join(args.chunks_dir, f"queries_{index}.fasta")).write_text(
            "".join(buffer), encoding="utf-8"
        )

    for _, record_id, sequence in queries:
        buffer.append(f">{record_id}\n{sequence}\n")
        if limit is not None and len(buffer) == limit:
            flush()
            index += 1
            buffer = []
    if buffer:
        flush()


def run_single(args):
    sequence = sample_queries(args.input)
    Path(args.output).write_text(f">{args.id}\n{sequence}\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Reduce sample inputs to colon-joined ColabFold query FASTAs"
    )
    sub = parser.add_subparsers(dest="command")

    single = sub.add_parser("single", help="reduce one sample input to its query FASTA")
    single.add_argument("input")
    single.add_argument("--id", required=True)
    single.add_argument("--output", required=True)

    batch = sub.add_parser("batch", help="convert and pre-batch a whole run")
    batch.add_argument("--manifest", required=True, help="TSV of sample id, record header, input file")
    batch.add_argument("--batch-size", type=int, default=-1, help="records per chunk FASTA (-1 = one chunk)")
    batch.add_argument("--chunks-dir", help="write pre-batched chunk FASTAs here")
    batch.add_argument("--queries-dir", help="write one query FASTA per sample here")

    args = parser.parse_args()
    if args.command is None:
        parser.print_help(file=sys.stderr)
        sys.exit(2)
    if args.command == "single":
        run_single(args)
    elif args.command == "batch":
        run_batch(args)


if __name__ == "__main__":
    main()
