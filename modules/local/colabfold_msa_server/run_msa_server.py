import argparse
import csv
from importlib.metadata import version
from pathlib import Path

from colabfold.colabfold import run_mmseqs2
from colabfold.input import msa_to_str, safe_filename


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--host-url", default="https://api.colabfold.com")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.input_csv.open(newline="") as csv_file:
        queries = list(csv.DictReader(csv_file))

    if not queries or not {"id", "sequence"}.issubset(queries[0]):
        raise ValueError("Input CSV must contain id and sequence columns")

    host_url = args.host_url.rstrip("/")
    user_agent = f"colabfold/{version('colabfold')}"
    for query in queries:
        sequence = query["sequence"].strip()
        query_sequences = [chain.strip() for chain in sequence.split(":")]
        if not sequence or any(not chain for chain in query_sequences):
            raise ValueError(f"Query {query['id']} has an empty sequence or chain")

        query_seqs_unique = []
        for chain in query_sequences:
            if chain not in query_seqs_unique:
                query_seqs_unique.append(chain)
        query_seqs_cardinality = [
            query_sequences.count(chain) for chain in query_seqs_unique
        ]
        prefix = args.output_dir / "msa" / safe_filename(query["id"])
        prefix.parent.mkdir(parents=True, exist_ok=True)

        unpaired_msa = run_mmseqs2(
            query_seqs_unique,
            str(prefix),
            use_env=True,
            host_url=host_url,
            user_agent=user_agent,
        )

        if len(query_seqs_unique) > 1:
            paired_msa = run_mmseqs2(
                query_seqs_unique,
                str(prefix),
                use_env=False,
                use_pairing=True,
                pairing_strategy="greedy",
                host_url=host_url,
                user_agent=user_agent,
            )
        elif query_seqs_cardinality[0] > 1:
            paired_msa = [
                f">{101 + index}\n{query_seqs_unique[0]}\n"
                for index in range(query_seqs_cardinality[0])
            ]
        else:
            paired_msa = None

        msa = msa_to_str(
            unpaired_msa,
            paired_msa,
            query_seqs_unique,
            query_seqs_cardinality,
        )
        (args.output_dir / f"{safe_filename(query['id'])}.a3m").write_text(msa)


if __name__ == "__main__":
    main()
