#!/usr/bin/env python3
"""Packed metric store for proteinfold.

Owns the on-disk format of the packed per-batch metric intermediates so that
the producer (pack_metrics.py) and the consumers (generate_multiqc_contents.py,
generate_report.py, populate_modelcif.py) never duplicate format logic.

Two files per inference batch replace the ~9 per-sample TSVs that
extract_metrics.py used to write:

``colabfold2_metrics.jsonl.gz``
    Scalar metrics, one JSON object per line::

        {"sample_id": "...", "rank": 0, "metric": "ptm",
         "chain_pair": null, "value": 0.841}

    ``metric`` is one of ``ptm``, ``iptm``, ``ipsae``, ``msa_depth``,
    ``msa_length``, ``chainwise_iptm``, ``chainwise_ipsae`` (the chainwise
    metrics carry a ``chain_pair`` label such as ``A:B``), plus the boolean
    flag ``plddt_structure_only``.

``colabfold2_metrics.npz``
    Per-residue / per-pair arrays, one compressed member per
    ``(sample_id, rank, kind)`` keyed ``"<sample_id>|<rank>|<kind>"``:

    * ``plddt``: 1-D float64 array of per-residue pLDDT.
    * ``pae``:   2-D float64 (n x n) PAE matrix.

    Values are quantized exactly like the legacy TSV formatting (pae to 4
    decimals) so consumers render values identically to the legacy TSVs.

The format is deliberately isolated here: switching to parquet (pyarrow) or
any other container later only means changing this module.
"""

import gzip
import io
import json
import os
import sys
import zipfile

import numpy as np

SCALARS_SUFFIX = "colabfold2_metrics.jsonl.gz"
ARRAYS_SUFFIX = "colabfold2_metrics.npz"



def array_key(sample_id, rank, kind):
    """npz member key for one (sample, rank, kind) array."""
    return f"{sample_id}|{rank}|{kind}"


def parse_array_key(key):
    """Inverse of :func:`array_key`; sample ids cannot contain '|'."""
    sample_id, rank, kind = key.split("|")
    return sample_id, int(rank), kind


def write_scalars(path, records):
    """Write scalar metric records to a gzipped JSONL file."""
    with gzip.open(path, "wt") as handle:
        for record in records:
            handle.write(json.dumps(record, allow_nan=False) + "\n")


def read_scalars(path):
    """Iterate scalar metric records from a gzipped JSONL file."""
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line:
                yield json.loads(line)


class ArrayStoreWriter:
    """Streaming writer for the compressed npz array store.

    Writes one .npy member per key into a deflated zip incrementally, so a
    whole-genome batch (``--inference_batch_size -1``) never holds more than
    one sample's arrays in memory. Members are read back transparently by
    :func:`read_arrays` (same layout ``np.savez_compressed`` produces).
    """

    def __init__(self, path):
        self._zip = zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED)

    def add(self, key, array):
        buffer = io.BytesIO()
        np.lib.format.write_array(buffer, np.asanyarray(array), allow_pickle=False)
        # Fixed member timestamp so identical content produces identical
        # bytes (zip default stamps the wall clock, which would make the
        # stores unstable for snapshot md5s and cache keys).
        info = zipfile.ZipInfo(f"{key}.npy", date_time=(1980, 1, 1, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        self._zip.writestr(info, buffer.getvalue())

    def close(self):
        self._zip.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
        return False


def write_arrays(path, arrays):
    """Write {key: array} to a compressed npz file in one shot."""
    with ArrayStoreWriter(path) as writer:
        for key, array in arrays.items():
            writer.add(key, array)


class ArrayStore:
    """Lazy read access to the compressed npz array store.

    Members are only decompressed when requested, so a whole-genome batch
    never materialises more than the arrays actually consumed.
    """

    def __init__(self, path):
        self._archive = np.load(path, allow_pickle=False)

    def get(self, sample_id, rank, kind):
        """Return the array for (sample_id, rank, kind), or None if absent."""
        key = array_key(sample_id, rank, kind)
        return self._archive[key] if key in self._archive.files else None

    def has(self, sample_id, kind):
        """True if any rank of `kind` exists for the sample."""
        return any(
            parse_array_key(key)[:1] == (sample_id,) and parse_array_key(key)[2] == kind
            for key in self._archive.files
        )

    def ranks(self, sample_id, kind):
        """Sorted ranks at which `kind` is stored for the sample."""
        found = {
            parse_array_key(key)[1]
            for key in self._archive.files
            if parse_array_key(key)[0] == sample_id and parse_array_key(key)[2] == kind
        }
        return sorted(found)

    def close(self):
        self._archive.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
        return False


def read_arrays(path):
    """Return {key: np.ndarray} from a compressed npz file (all members)."""
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def quantize_pae(matrix):
    """Quantize a PAE matrix exactly like the legacy ``{:.4f}`` TSV cells."""
    return np.vectorize(lambda v: float(f"{v:.4f}"), otypes=[np.float64])(np.asarray(matrix, dtype=np.float64))


def read_msa_depth(a3m_path):
    """MSA depth, counted exactly like extract_metrics.a3m_to_int does.

    Counts sequence lines (non-``>`` lines), skipping a leading non-FASTA
    cardinality line such as ColabFold's ``#N<TAB>M``. Avoids materialising
    the whole alignment, which matters for deep MSAs at genome scale.
    """
    depth = 0
    first = True
    with open(a3m_path) as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if first:
                first = False
                if stripped and not stripped.startswith(">"):
                    continue
            if not stripped.startswith(">"):
                if stripped:
                    depth += 1
    return depth


def read_msa_length(a3m_path):
    """Alignment length: insert-stripped length of the first sequence."""
    first = True
    with open(a3m_path) as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if first:
                first = False
                if stripped and not stripped.startswith(">"):
                    continue
            if stripped.startswith(">"):
                continue
            if stripped:
                return len("".join(char for char in stripped if not char.islower()))
    return 0
