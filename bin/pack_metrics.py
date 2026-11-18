#!/usr/bin/env python3
"""Pack per-sample inference metrics into per-batch packed intermediates.

Replaces the per-sample TSVs that extract_metrics.py writes (ptm/iptm/ipsae,
chainwise scores, per-rank PAE matrices, per-residue pLDDT, integer-encoded
MSA) with two files per batch:

* ``colabfold2_metrics.jsonl.gz`` -- scalar metrics (see metrics_store.py)
* ``colabfold2_metrics.npz``     -- pLDDT vectors and PAE matrices

Values are quantized like the legacy TSV formatting so consumers render
them identically. The full MSA is NOT packed: consumers read the a3m
directly (the integer-encoded TSV was a size-multiplied duplicate).

Runs once per inference batch task over all samples in ``--ids``:

    pack_metrics.py --ids sample_ids.txt --struct-ext cif --raw-dir raw \\
        --out-scalars colabfold2_metrics.jsonl.gz \\
        --out-arrays colabfold2_metrics.npz
"""

import argparse
import glob
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import extract_metrics  # noqa: E402  (staged alongside by Nextflow)
import metrics_store  # noqa: E402


def score_files_for(raw_dir, sample_id):
    return sorted(glob.glob(os.path.join(raw_dir, f"{sample_id}_scores_rank_*.json")))


def struct_files_for(raw_dir, sample_id, struct_ext):
    return sorted(glob.glob(os.path.join(raw_dir, f"{sample_id}_rank_*.{struct_ext}")))


def pack_sample(sample_id, raw_dir, struct_ext, scalars_out, arrays_out):
    """Pack one sample's metrics into the open store writers."""
    score_files = score_files_for(raw_dir, sample_id)
    structs = struct_files_for(raw_dir, sample_id, struct_ext)
    if not score_files:
        print(f"pack_metrics: no scores JSON for {sample_id}; skipping metric extraction", file=sys.stderr)
        return

    records = []
    collected = extract_metrics.collect_colabfold_metrics(score_files, name=sample_id)

    for rank_str, value in collected["ptm_rows"]:
        records.append({"sample_id": sample_id, "rank": int(rank_str), "metric": "ptm", "chain_pair": None, "value": value})
    # Legacy derived iptm values are pre-formatted strings; keep them as
    # strings so rehydration renders them verbatim.
    for rank_str, value in collected["iptm_rows"]:
        records.append(
            {"sample_id": sample_id, "rank": int(rank_str), "metric": "iptm", "chain_pair": None, "value": value}
        )
    for rank_str, value in collected["ipsae_rows"]:
        records.append(
            {"sample_id": sample_id, "rank": int(rank_str), "metric": "ipsae", "chain_pair": None, "value": value}
        )
    for metric, chainwise in (
        ("chainwise_iptm", collected["chainwise_iptm"]),
        ("chainwise_ipsae", collected["chainwise_ipsae"]),
    ):
        for rank, pair_scores in chainwise.items():
            for pair_label, value in pair_scores:
                records.append(
                    {
                        "sample_id": sample_id,
                        "rank": int(rank),
                        "metric": metric,
                        "chain_pair": pair_label,
                        "value": value,
                    }
                )

    if structs:
        computed = extract_metrics.compute_plddt_columns(structs, score_files, name=sample_id)
        if computed["structure_only"]:
            records.append(
                {"sample_id": sample_id, "rank": 0, "metric": "plddt_structure_only", "chain_pair": None, "value": 1}
            )
        # Store the columns even when structure-only: rehydration writes the
        # legacy 'n/a' table from the array length, matching extract_metrics.
        for rank, column in zip(computed["ranks"], computed["plddt_cols"]):
            arrays_out.add(metrics_store.array_key(sample_id, rank, "plddt"), np.asarray(column, dtype=np.float64))

    for rank, pae in collected["pae_by_rank"].items():
        arrays_out.add(metrics_store.array_key(sample_id, rank, "pae"), metrics_store.quantize_pae(pae))

    a3m_path = os.path.join(raw_dir, f"{sample_id}.a3m")
    if os.path.exists(a3m_path):
        records.append(
            {"sample_id": sample_id, "rank": 0, "metric": "msa_depth", "chain_pair": None, "value": metrics_store.read_msa_depth(a3m_path)}
        )
        records.append(
            {"sample_id": sample_id, "rank": 0, "metric": "msa_length", "chain_pair": None, "value": metrics_store.read_msa_length(a3m_path)}
        )

    for record in records:
        scalars_out.write(json.dumps(record, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", required=True, help="File with one sanitized sample id per line")
    parser.add_argument("--struct-ext", default="cif", choices=["cif", "pdb"], help="Structure file extension to glob")
    parser.add_argument("--raw-dir", default="raw", help="Directory holding colabfold_batch outputs")
    parser.add_argument("--out-scalars", default=metrics_store.SCALARS_SUFFIX)
    parser.add_argument("--out-arrays", default=metrics_store.ARRAYS_SUFFIX)
    args = parser.parse_args()

    with open(args.ids) as handle:
        sample_ids = [line.strip() for line in handle if line.strip()]

    # Stream: one sample's arrays in memory at a time, so a whole-genome
    # single-batch run (--inference_batch_size -1) stays bounded.
    with gzip.open(args.out_scalars, "wt") as scalars_out, metrics_store.ArrayStoreWriter(args.out_arrays) as arrays_out:
        for sample_id in sample_ids:
            pack_sample(sample_id, args.raw_dir, args.struct_ext, scalars_out, arrays_out)

    print(f"pack_metrics: packed {len(sample_ids)} sample(s) into {args.out_scalars} and {args.out_arrays}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
