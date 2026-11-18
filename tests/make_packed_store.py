#!/usr/bin/env python3
"""Build a packed metric store pair for the assemble_modelcif tests.

Reads the legacy per-sample metric TSVs (the same fixtures the per-sample
MODELCIF_ASSEMBLE tests used) and packs them into the
``<prefix>_metrics.jsonl.gz`` / ``<prefix>_metrics.npz`` pair consumed by
``populate_modelcif.py --manifest``.

Usage:
    make_packed_store.py --sample-id test --out-prefix path/prefix \
        --plddt rank_0=plddt0.tsv,rank_1=plddt1.tsv \
        --pae rank_0=pae0.tsv,rank_1=pae1.tsv \
        --scores rank_0=0.5,rank_1=0.6 \
        [--chainwise-iptm rank_0=A:B=0.1,...]
"""

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "bin"))

import metrics_store  # noqa: E402


def parse_rank_spec(spec):
    """Parse 'rank_0=a,rank_1=b' into an ordered {rank: value} dict."""
    if not spec:
        return {}
    out = {}
    for item in spec.split(","):
        rank, _, value = item.partition("=")
        out[int(rank.split("_", 1)[1])] = value
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--plddt", default="", help="rank_<n>=<tsv path with one float per residue>")
    parser.add_argument("--pae", default="", help="rank_<n>=<tsv matrix path>")
    parser.add_argument("--scores", default="", help="rank_<n>=<float> shared by ptm/iptm")
    parser.add_argument("--chainwise-iptm", default="", help="rank_<n>=<chainA>:<chainB>=<float>,...")
    parser.add_argument("--chainwise-ipsae", default="", help="rank_<n>=<chainA>:<chainB>=<float>,...")
    args = parser.parse_args()

    plddt = parse_rank_spec(args.plddt)
    pae = parse_rank_spec(args.pae)
    scores = parse_rank_spec(args.scores)
    chainwise_iptm = parse_rank_spec(args.chainwise_iptm)
    chainwise_ipsae = parse_rank_spec(args.chainwise_ipsae)

    scalars_path = f"{args.out_prefix}_metrics.jsonl.gz"
    arrays_path = f"{args.out_prefix}_metrics.npz"

    records = []
    for rank, value in sorted(scores.items()):
        records.append({
            "sample_id": args.sample_id, "rank": rank, "metric": "ptm",
            "chain_pair": None, "value": float(value),
        })
        records.append({
            "sample_id": args.sample_id, "rank": rank, "metric": "iptm",
            "chain_pair": None, "value": float(value),
        })
    for rank, spec in sorted(chainwise_iptm.items()):
        for pair, value in (item.split("=") for item in spec.split(";")):
            records.append({
                "sample_id": args.sample_id, "rank": rank, "metric": "chainwise_iptm",
                "chain_pair": pair, "value": float(value),
            })
    for rank, spec in sorted(chainwise_ipsae.items()):
        for pair, value in (item.split("=") for item in spec.split(";")):
            records.append({
                "sample_id": args.sample_id, "rank": rank, "metric": "chainwise_ipsae",
                "chain_pair": pair, "value": float(value),
            })
    metrics_store.write_scalars(scalars_path, records)

    arrays = {}
    for rank, path in sorted(plddt.items()):
        values = [float(line) for line in Path(path).read_text().split() if line.strip()]
        arrays[metrics_store.array_key(args.sample_id, rank, "plddt")] = values
    for rank, path in sorted(pae.items()):
        matrix = [
            [float(cell) for cell in line.split()] if "\t" in line else [float(line)]
            for line in Path(path).read_text().splitlines() if line.strip()
        ]
        arrays[metrics_store.array_key(args.sample_id, rank, "pae")] = matrix
    metrics_store.write_arrays(arrays_path, arrays)

    print(scalars_path, arrays_path)


if __name__ == "__main__":
    main()
