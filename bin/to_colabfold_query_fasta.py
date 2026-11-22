#!/usr/bin/env python3
"""Reduce sample inputs to colon-joined ColabFold query FASTAs and run the
MSA cache lookups/stores around them.

colabfold_search treats a sequence containing ':' as a complex, so every
chain of a sample is colon-joined into one query record: multi-record
FASTAs are joined across their records, and AlphaFold3 JSONs are reduced
to their protein chains (repeated per copy, in entity order), so the
search sees the same chains the JSON describes. Non-protein entities
(RNA, DNA, ligands) are ignored here; for JSON inputs the original
document is repopulated with the generated MSAs downstream
(merge_alphafold3_msa.py).

The MSA cache lives here too, because both sides of a cached search must
hash the query content identically: the conversion looks queries up
before the search, and the search stores its results under the same
content keys. Every entry is a zip archive holding the AF3-format MSA
JSON and its a3m, keyed by the SHA-256 of the query sequences plus the
database, search tool and search arguments.

Subcommands:
  single INPUT --id ID --output PATH
      Reduce one sample input to its single-query FASTA.
  batch --manifest TSV [--batch-size N] (--chunks-dir D | --queries-dir D)
        [--hits-dir D] [--cache-dir D --cache-id X --tool-id Y --cache-args A]
      Process a whole run in one go without per-sample intermediates: each
      record is cache-looked-up in memory, hits are extracted straight to
      the cache-result directory, and the misses are concatenated on the
      fly into chunk FASTAs of --batch-size records (or written one per
      sample with --queries-dir).
  cache-lookup QUERY <cache_id> <tool_id> <args> <dir> <destination> <artifact>
      Extract the cached artifact (af3_json or a3m) for the query into
      <destination>; writes nothing on a miss -- the written file is the
      hit signal, not the exit code.
  cache-store QUERY <cache_id> <tool_id> <args> <dir> <json> <a3m>
      Store one query's results under its content key.
  cache-store-chunk CHUNK <cache_id> <tool_id> <args> <dir> <results_dir>
      Store every record of a concatenated query chunk separately: the
      record header names the result files (sanitised the same way
      colabfold_search names them), and each record is keyed by its own
      sequence, so entries match the per-sample lookups exactly.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
import zipfile
import zlib
from pathlib import Path

CACHE_FORMAT_VERSION = 3

# colabfold_search sanitises record headers into result file names with the
# same character class; keep both sides identical.
_SANITIZE_RE = re.compile(r"[^A-Za-z0-9_.-]")


def sanitise_name(name):
    return _SANITIZE_RE.sub("_", name)


# --------------------------------------------------------------------------
# MSA cache
# --------------------------------------------------------------------------


def key_for_sequences(sequences, cache_id, tool_id, args=""):
    payload = {
        "sequences": sequences,
        "database": cache_id or "",
        "search_tool": tool_id or "",
        "args": args,
        "format": CACHE_FORMAT_VERSION,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def parse_query_file(path):
    """Read a query FASTA or id,sequence CSV as a list of sequences."""
    sequences = []
    with open(path, encoding="utf-8", newline="") as handle:
        if str(path).lower().endswith(".csv"):
            reader = csv.DictReader(handle)
            if reader.fieldnames != ["id", "sequence"]:
                raise ValueError("Expected id,sequence query CSV")
            sequences = [row["sequence"] for row in reader]
        else:
            for line in handle:
                if line.startswith(">"):
                    sequences.append("")
                elif sequences:
                    sequences[-1] += line.strip()
    return sequences


def key_for(query, cache_id, tool_id, args=""):
    return key_for_sequences(parse_query_file(query), cache_id, tool_id, args)


def valid_json(path):
    try:
        with open(path, encoding="utf-8") as handle:
            value = json.load(handle)
        return (
            isinstance(value, dict)
            and isinstance(value.get("sequences"), list)
            and bool(value["sequences"])
            and all(isinstance(item, dict) for item in value["sequences"])
            and all(any(kind in item for kind in ("protein", "rna", "dna", "ligand")) for item in value["sequences"])
            and all(
                isinstance(protein, dict)
                and isinstance(protein.get("id"), (str, list))
                and isinstance(protein.get("sequence"), str)
                and bool(protein["sequence"])
                and isinstance(protein.get("unpairedMsa"), str)
                and isinstance(protein.get("pairedMsa"), str)
                for item in value["sequences"]
                for protein in ([item["protein"]] if "protein" in item else [])
            )
        )
    except (OSError, ValueError, TypeError):
        return False


def valid_a3m(path):
    try:
        with open(path, encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]
        while lines and lines[0].startswith("#"):
            if len(lines[0]) == 1:
                return False
            lines.pop(0)
        if len(lines) < 2:
            return False
        has_header = False
        has_sequence = False
        for line in lines:
            if line.startswith(">"):
                if (has_header and not has_sequence) or len(line) == 1:
                    return False
                has_header = True
                has_sequence = False
            elif not has_header or line.startswith("#") or not all(c.isalpha() or c in "-.*:" for c in line):
                return False
            else:
                has_sequence = True
        return has_header and has_sequence
    except (OSError, UnicodeError):
        return False


def cache_path(directory, key):
    if not os.path.isabs(directory):
        raise ValueError("MSA cache directory must be an absolute path")
    return os.path.join(directory, key + ".zip")


def lookup(directory, key, destination, artifact):
    members = {"af3_json": ("alignment.json", valid_json), "a3m": ("alignment.a3m", valid_a3m)}
    member, validate = members[artifact]
    source = cache_path(directory, key)
    temporary = None
    try:
        with zipfile.ZipFile(source) as archive:
            fd, temporary = tempfile.mkstemp(prefix=".msa-hit-", dir=os.path.dirname(os.path.abspath(destination)))
            with os.fdopen(fd, "wb") as output, archive.open(member) as content:
                while True:
                    chunk = content.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
            if not validate(temporary):
                return False
            os.replace(temporary, destination)
            return True
    except (OSError, ValueError, KeyError, EOFError, RuntimeError, zipfile.BadZipFile, zlib.error):
        return False
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)


def store(directory, key, json_source, a3m_source):
    if not valid_json(json_source) or not valid_a3m(a3m_source):
        raise ValueError("MMSEQS output is not a valid MSA artifact")
    destination = cache_path(directory, key)
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=key + ".", suffix=".tmp", dir=directory)
    try:
        with os.fdopen(fd, "wb") as output:
            with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.write(a3m_source, "alignment.a3m")
                archive.write(json_source, "alignment.json")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def fasta_records(path):
    """Read a (concatenated) FASTA as one (header, sequence) pair per record."""
    records = []
    header = None
    current = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(current)))
                header = line[1:].strip()
                current = []
            elif line:
                if header is None:
                    raise ValueError(f"FASTA input {path} has sequence data before any header line")
                current.append(line)
    if header is not None:
        records.append((header, "".join(current)))
    return records


def store_chunk(query, cache_id, tool_id, args, directory, results_dir):
    for header, sequence in fasta_records(query):
        sample_id = sanitise_name(header)
        json_source = os.path.join(results_dir, sample_id + ".json")
        a3m_source = os.path.join(results_dir, sample_id + ".a3m")
        key = key_for_sequences([sequence], cache_id, tool_id, args)
        try:
            store(directory, key, json_source, a3m_source)
        except (OSError, ValueError) as error:
            print(f"Warning: could not store MSA cache entry for {sample_id}: {error}", file=sys.stderr)


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
    cache_enabled = args.cache_dir is not None
    if cache_enabled and not os.path.isabs(args.cache_dir):
        raise ValueError("MSA cache directory must be an absolute path")

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

    if cache_enabled:
        os.makedirs(args.hits_dir, exist_ok=True)

    misses = []
    for sample_id, record_id, sequence in entries:
        if cache_enabled:
            key = key_for_sequences([sequence], args.cache_id, args.tool_id, args.cache_args)
            destination = os.path.join(args.hits_dir, sample_id + ".json")
            if lookup(args.cache_dir, key, destination, "af3_json"):
                # The cached AF3-format MSA JSON stands in for the search.
                continue
        misses.append((sample_id, record_id, sequence))

    if args.queries_dir:
        os.makedirs(args.queries_dir, exist_ok=True)
        for sample_id, record_id, sequence in misses:
            Path(os.path.join(args.queries_dir, sample_id + ".fasta")).write_text(
                f">{record_id}\n{sequence}\n", encoding="utf-8"
            )
        return

    os.makedirs(args.chunks_dir, exist_ok=True)
    limit = args.batch_size if args.batch_size > 0 else None
    index = 0
    buffer = []

    def flush():
        Path(os.path.join(args.chunks_dir, f"queries_{index}.fasta")).write_text(
            "".join(buffer), encoding="utf-8"
        )

    for _, record_id, sequence in misses:
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
        description="Reduce sample inputs to ColabFold query FASTAs and run the MSA cache around them"
    )
    sub = parser.add_subparsers(dest="command")

    single = sub.add_parser("single", help="reduce one sample input to its query FASTA")
    single.add_argument("input")
    single.add_argument("--id", required=True)
    single.add_argument("--output", required=True)

    batch = sub.add_parser("batch", help="convert, cache-look-up and pre-batch a whole run")
    batch.add_argument("--manifest", required=True, help="TSV of sample id, record header, input file")
    batch.add_argument("--batch-size", type=int, default=-1, help="records per chunk FASTA (-1 = one chunk)")
    batch.add_argument("--chunks-dir", help="write pre-batched chunk FASTAs here")
    batch.add_argument("--queries-dir", help="write one query FASTA per sample here")
    batch.add_argument("--hits-dir", default="cache_result", help="directory for cached MSA JSONs")
    batch.add_argument("--cache-dir", help="MSA cache directory (omit to disable lookups)")
    batch.add_argument("--cache-id", default="", help="cache database identifier")
    batch.add_argument("--tool-id", default="", help="cache search tool identifier")
    batch.add_argument("--cache-args", default="", help="search arguments hashed into the key")

    lookup_cmd = sub.add_parser("cache-lookup", help="extract a cached MSA artifact for a query")
    lookup_cmd.add_argument("query")
    lookup_cmd.add_argument("cache_id")
    lookup_cmd.add_argument("tool_id")
    lookup_cmd.add_argument("args")
    lookup_cmd.add_argument("directory")
    lookup_cmd.add_argument("destination")
    lookup_cmd.add_argument("artifact", choices=["af3_json", "a3m"])

    store_cmd = sub.add_parser("cache-store", help="store one query's MSA results")
    store_cmd.add_argument("query")
    store_cmd.add_argument("cache_id")
    store_cmd.add_argument("tool_id")
    store_cmd.add_argument("args")
    store_cmd.add_argument("directory")
    store_cmd.add_argument("json_source")
    store_cmd.add_argument("a3m_source")

    store_chunk_cmd = sub.add_parser("cache-store-chunk", help="store every record of a concatenated query")
    store_chunk_cmd.add_argument("query")
    store_chunk_cmd.add_argument("cache_id")
    store_chunk_cmd.add_argument("tool_id")
    store_chunk_cmd.add_argument("args")
    store_chunk_cmd.add_argument("directory")
    store_chunk_cmd.add_argument("results_dir")

    args = parser.parse_args()
    if args.command is None:
        parser.print_help(file=sys.stderr)
        sys.exit(2)
    if args.command == "single":
        run_single(args)
    elif args.command == "batch":
        run_batch(args)
    elif args.command == "cache-lookup":
        lookup(args.directory, key_for(args.query, args.cache_id, args.tool_id, args.args), args.destination, args.artifact)
    elif args.command == "cache-store":
        try:
            store(args.directory, key_for(args.query, args.cache_id, args.tool_id, args.args), args.json_source, args.a3m_source)
        except (OSError, ValueError) as error:
            print(f"Warning: could not store MSA cache: {error}", file=sys.stderr)
    elif args.command == "cache-store-chunk":
        store_chunk(args.query, args.cache_id, args.tool_id, args.args, args.directory, args.results_dir)


if __name__ == "__main__":
    main()
