#!/usr/bin/env python3
"""Generate MultiQC custom content for one prediction mode.

Reads the routed metric TSVs and writes proteinfold_<mode>_generalstats_mqc.json
and proteinfold_<mode>_plddt_lineplot_mqc.json for the stock MULTIQC module.
"""

import argparse
import array
import ast
import csv
import gzip
import json
import math
import re
import struct
import sys
import zipfile
from pathlib import Path

MODE_LABELS = {
    "alphafold2": "AlphaFold2",
    "alphafold3": "AlphaFold3",
    "colabfold": "ColabFold",
    "esmfold": "ESMFold",
    "boltz": "Boltz",
}

RANK_COLUMN_RE = re.compile(r"^rank_(\d+)$")




def warn(message):
    """Print a warning to stderr; warnings never abort the report."""
    print(f"WARNING: {message}", file=sys.stderr)


def fail(message):
    """Print an error and exit non-zero (pipeline wiring problems only)."""
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(2)


def clean_float(value):
    """Return a finite float, or None for missing/NaN/infinite values."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def set_metric(rows, sample_name, column, value, source):
    """Set a metric, warning on duplicate-id collapse (the samplesheet does not enforce unique ids)."""
    if column in rows.setdefault(sample_name, {}):
        warn(
            f"Duplicate sample data: '{sample_name}' already has '{column}' values; "
            f"values from {Path(source).name} overwrite the earlier file. "
            "The input samplesheet schema does not enforce unique ids."
        )
    rows[sample_name][column] = value


# ---------------------------------------------------------------------------
# Packed metric stores (bin/pack_metrics.py output)
# ---------------------------------------------------------------------------

SCALARS_MARKER = "_colabfold2_metrics.jsonl.gz"
ARRAYS_MARKER = "_colabfold2_metrics.npz"


def is_packed_scalars(name):
    return name.endswith(".jsonl.gz") and "colabfold2_metrics" in name


def is_packed_arrays(name):
    return name.endswith(".npz") and "colabfold2_metrics" in name


def read_npz_float64(path):
    """Minimal stdlib reader for the packed npz store (float64 members only).

    The packed arrays are plain ``.npy`` v1.0 members inside a zip; this
    avoids a numpy dependency in the MultiQC content generator container.
    Returns ``{key: list}`` (1-D) or ``{key: list[list]}`` (2-D).
    """
    unpacked = {}
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if not name.endswith(".npy"):
                continue
            raw = archive.read(name)
            if raw[:6] != b"\x93NUMPY":
                raise ValueError(f"{path.name}: member {name} is not an npy file")
            major = raw[6]
            if major != 1:
                raise ValueError(f"{path.name}: member {name}: unsupported npy version {major}")
            (header_len,) = struct.unpack("<H", raw[8:10])
            header = ast.literal_eval(raw[10 : 10 + header_len].decode("latin1"))
            if header["descr"] != "<f8":
                raise ValueError(f"{path.name}: member {name}: unsupported dtype {header['descr']}")
            shape = header["shape"]
            payload = raw[10 + header_len :]
            values = array.array("d")
            values.frombytes(payload)
            if sys.byteorder != "little":
                values.byteswap()
            flat = values.tolist()
            if len(shape) == 1:
                unpacked[name[: -len(".npy")]] = flat
            elif len(shape) == 2:
                rows_n, cols_n = shape
                unpacked[name[: -len(".npy")]] = [flat[i * cols_n : (i + 1) * cols_n] for i in range(rows_n)]
            else:
                raise ValueError(f"{path.name}: member {name}: unsupported shape {shape}")
    return unpacked


def parse_packed_pair(scalars_path, arrays_path, rows, line_series, mode_label):
    """Fold one batch's packed metric stores into the generalstats rows."""
    parsed_records = 0
    rank_scores = {}       # sample_id -> metric -> {rank: value}
    structure_only = set()
    with gzip.open(scalars_path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            parsed_records += 1
            metric = record["metric"]
            sample_id = record["sample_id"]
            sample_name = f"{sample_id}_{mode_label}"
            if metric == "plddt_structure_only":
                if record["value"]:
                    structure_only.add(sample_id)
            elif metric == "msa_depth":
                if record["value"]:
                    set_metric(rows, sample_name, "msa_depth", record["value"], scalars_path)
            elif metric in ("ptm", "iptm"):
                rank_scores.setdefault(sample_id, {}).setdefault(metric, {})[record["rank"]] = record["value"]

    for sample_id, metrics in rank_scores.items():
        sample_name = f"{sample_id}_{mode_label}"
        for metric, scores in metrics.items():
            cleaned = {rank: value for rank, value in scores.items() if clean_float(value) is not None}
            if not cleaned:
                warn(f"{scalars_path.name}: {metric} for {sample_id} yielded no usable values; skipping")
                continue
            top_rank = min(cleaned)
            set_metric(rows, sample_name, metric, cleaned[top_rank], scalars_path)
            for rank, value in cleaned.items():
                if rank != top_rank:
                    set_metric(rows, f"{sample_name}_rank_{rank}", metric, value, scalars_path)

    arrays = read_npz_float64(arrays_path)
    plddt_by_sample = {}   # sample_id -> {rank: values}
    for key, values in arrays.items():
        sample_id, rank_str, kind = key.rsplit("|", 2)
        if kind != "plddt":
            continue  # PAE matrices stay in the detailed report viewer, like legacy
        if sample_id in structure_only:
            continue
        plddt_by_sample.setdefault(sample_id, {})[int(rank_str)] = values

    for sample_id, rank_values in plddt_by_sample.items():
        sample_name = f"{sample_id}_{mode_label}"
        rank_means = {}
        rank_points = {}
        for rank, values in rank_values.items():
            # The packed arrays are plain per-residue lists whose first entry
            # is residue 1; the legacy TSV path carried an explicit 1-based
            # Positions column, so keep the same numbering here.
            points = [[position, value] for position, value in enumerate(values, start=1) if math.isfinite(value)]
            if points:
                rank_points[rank] = points
                rank_means[rank] = sum(value for _, value in points) / len(points)
        if not rank_means:
            continue
        top_rank = min(rank_means)
        set_metric(rows, sample_name, "mean_plddt", rank_means[top_rank], arrays_path)
        for rank, mean in rank_means.items():
            if rank != top_rank:
                set_metric(rows, f"{sample_name}_rank_{rank}", "mean_plddt", mean, arrays_path)
        if rank_points:
            line_series.setdefault(sample_name, {})
            for rank, points in rank_points.items():
                line_series[sample_name][f"rank_{rank}"] = points

    return parsed_records


GENERALSTATS_HEADERS = {
    "msa_depth": {
        "title": "Related sequence depth (MSA)",
        "description": "The number of related sequences (across the whole protein) that could be retrieved from the "
        "MSA (Multiple Sequence Alignment) stage",
        "namespace": "proteinfold",
        "format": "{:,.0f}",
    },
    "mean_plddt": {
        "title": "Structure confidence (average pLDDT)",
        "description": "Structure prediction confidence score across all residues in the top ranked protein "
        "structure - from the mean pLDDT (predicted Local Distance Difference Test) value",
        "namespace": "proteinfold",
        "max": 100,
        "min": 0,
        "cond_formatting_rules": {
            "very-low": [{"lt": 50}],
            "low": [{"gt": 50}, {"lt": 70}],
            "high": [{"gt": 70}, {"lt": 90}],
            "very-high": [{"gt": 90}],
        },
        "cond_formatting_colours": [
            {"very-low": "#f0743e"},
            {"low": "#f9d613"},
            {"high": "#60c2e8"},
            {"very-high": "#014ecc"},
        ],
    },
    "iptm": {
        "title": "Interface accuracy (ipTM)",
        "description": "Accuracy of the relative positions of two protein subunits from a multimer calculation - "
        "from the ipTM (interface predicted Template Modelling) score",
        "namespace": "proteinfold",
        "max": 1,
        "min": 0,
        "format": "{:,.2f}",
        "scale": "Purples",
    },
    "ptm": {
        "title": "Global accuracy (TM)",
        "description": "Global accuracy of the protein folded, less sensitive to localised inaccuracies than raw 3D "
        "atomic deviations (RMSD) - from the pTM (predicted Template Modelling) score",
        "namespace": "proteinfold",
        "max": 1,
        "min": 0,
        "format": "{:,.2f}",
        "scale": "Blues",
    },
}

LINEGRAPH_PCONFIG = {
    "id": "proteinfold_plddt_lineplot",
    "title": "ProteinFold: pLDDT by Position",
    "xlab": "Residue Position",
    "ylab": "pLDDT Score",
    "ymin": 0,
    "ymax": 100,
}


def dump_json(path, payload):
    """Write JSON, refusing to emit NaN/Infinity tokens (allow_nan=False)."""
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False, sort_keys=False)
        handle.write("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate MultiQC custom-content JSON from proteinfold metric TSVs."
    )
    parser.add_argument(
        "--model",
        required=True,
        help="Prediction-mode key from the pipeline meta map (e.g. alphafold2, boltz); "
        "sample names are qualified with this mode, never with 'UNKNOWN'.",
    )
    parser.add_argument("--output-dir", default=".", help="Directory to write the *_mqc.json files into")
    parser.add_argument("metric_files", nargs="+", help="Packed metric stores (*_colabfold2_metrics.jsonl.gz / .npz)")
    args = parser.parse_args(argv)

    model_key = (args.model or "").strip()
    if not model_key:
        fail("No prediction mode given (--model); refusing to emit '_UNKNOWN' sample names")
    mode_label = MODE_LABELS.get(model_key, model_key)

    rows = {}
    line_series = {}
    parsed_files = 0

    packed_scalars = [p for p in args.metric_files if is_packed_scalars(Path(p).name)]
    packed_arrays = [p for p in args.metric_files if is_packed_arrays(Path(p).name)]
    if not packed_scalars and not packed_arrays:
        fail("No packed metric stores were passed; the MULTIQC input channel is miswired")
    arrays_by_base = {Path(p).name[: -len(ARRAYS_MARKER)]: p for p in packed_arrays}
    for raw_path in packed_scalars:
        filepath = Path(raw_path)
        base = filepath.name[: -len(SCALARS_MARKER)]
        arrays_path = arrays_by_base.get(base)
        if arrays_path is None:
            fail(f"Packed scalar store {filepath.name} has no matching array store '{base}{ARRAYS_MARKER}'")
        parsed_files += parse_packed_pair(filepath, Path(arrays_path), rows, line_series, mode_label)

    metrics_seen = {column for data in rows.values() for column in data}
    headers = {column: dict(GENERALSTATS_HEADERS[column]) for column in GENERALSTATS_HEADERS if column in metrics_seen}
    generalstats = {
        "id": "proteinfold",
        "section_name": "ProteinFold",
        "description": "Summary metrics for protein structure prediction "
        "(average pLDDT, MSA depth, pTM and ipTM per prediction)",
        "plot_type": "generalstats",
        "headers": headers,
        "data": rows,
    }

    datasets = []
    data_labels = []
    for sample_name in sorted(line_series):
        if not line_series[sample_name]:
            continue
        datasets.append(line_series[sample_name])
        data_labels.append({"name": sample_name, "ylab": "pLDDT score"})

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    generalstats_path = output_dir / f"proteinfold_{model_key}_generalstats_mqc.json"
    dump_json(generalstats_path, generalstats)

    linegraph_path = None
    if datasets:
        pconfig = dict(LINEGRAPH_PCONFIG)
        pconfig["data_labels"] = data_labels
        linegraph = {
            "id": "proteinfold_plddt_lineplot",
            "section_name": "ProteinFold: pLDDT by Position",
            "description": "Per-residue confidence scores across all predicted ranks "
            "(switch between predictions with the tabs above the plot)",
            "plot_type": "linegraph",
            "pconfig": pconfig,
            "data": datasets,
        }
        linegraph_path = output_dir / f"proteinfold_{model_key}_plddt_lineplot_mqc.json"
        dump_json(linegraph_path, linegraph)

    if linegraph_path is None:
        print(
            f"Generated {generalstats_path.name} from {parsed_files} metric file(s): "
            f"{len(rows)} sample row(s), no per-residue pLDDT series"
        )
    else:
        print(
            f"Generated {generalstats_path.name} and {linegraph_path.name} "
            f"from {parsed_files} metric file(s): {len(rows)} sample row(s), "
            f"{len(datasets)} pLDDT dataset(s), "
            f"metrics seen: {', '.join(sorted(metrics_seen)) or 'none'}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
