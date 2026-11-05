#!/usr/bin/env python3
import csv
import hashlib
import json
import os
import tempfile
import sys
import zipfile

CACHE_FORMAT_VERSION = 3


def key_for(fasta, cache_id, tool_id, args=""):
    sequences = []
    with open(fasta, encoding="utf-8", newline="") as handle:
        if str(fasta).lower().endswith('.csv'):
            reader = csv.DictReader(handle)
            if reader.fieldnames != ['id', 'sequence']:
                raise ValueError("Expected id,sequence query CSV")
            sequences = [row['sequence'] for row in reader]
        else:
            for line in handle:
                if line.startswith(">"):
                    sequences.append("")
                elif sequences:
                    sequences[-1] += line.strip()
    payload = {"sequences": sequences, "database": cache_id or "", "search_tool": tool_id or "", "args": args, "format": CACHE_FORMAT_VERSION}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


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
        while lines and lines[0].startswith('#'):
            if len(lines[0]) == 1:
                return False
            lines.pop(0)
        if len(lines) < 2:
            return False
        has_header = False
        has_sequence = False
        for line in lines:
            if line.startswith('>'):
                if (has_header and not has_sequence) or len(line) == 1:
                    return False
                has_header = True
                has_sequence = False
            elif not has_header or line.startswith('#') or not all(c.isalpha() or c in '-.*:' for c in line):
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
                while chunk := content.read(1024 * 1024):
                    output.write(chunk)
            if not validate(temporary):
                return False
            os.replace(temporary, destination)
            return True
    except (OSError, ValueError, KeyError, EOFError, RuntimeError, zipfile.BadZipFile):
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


if __name__ == "__main__":
    action, fasta, cache_id, tool_id, args, directory, *paths = sys.argv[1:]
    key = key_for(fasta, cache_id, tool_id, args)
    if action == "lookup":
        destination, artifact = paths
        lookup(directory, key, destination, artifact)
    elif action == "store":
        json_source, a3m_source = paths
        try:
            store(directory, key, json_source, a3m_source)
        except (OSError, ValueError) as error:
            print(f"Warning: could not store MSA cache: {error}", file=sys.stderr)
    else:
        raise ValueError(action)
