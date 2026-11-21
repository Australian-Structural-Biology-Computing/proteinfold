#!/usr/bin/env python

import math
import os
import argparse
import csv
import json
from matplotlib import pyplot as plt
import numpy as np
import base64
import re
from Bio import PDB

def is_missing_input(path):
    """True when a metric input is absent or empty.

    The inference task signals an unavailable MSA (the ESMFold2
    language-model variants fold without one) with a zero-size TSV, which
    this treats the same as a missing file.
    """
    return not path or not os.path.exists(path) or os.path.getsize(path) == 0
# PAE values above this cap are clipped before quantisation; the report viewer
# uses the same cap, so displayed values and colours stay consistent.
PAE_MAX = 30.0
# Maximum number of MSA rows embedded in the interactive report. Rows are
# sorted by sequence identity before sampling, so the identity gradient of the
# coverage plot is preserved even for very deep MSAs.
MSA_MAX_ROWS = 1024

def encode_u8_b64(values):
    """Quantise floats in [0, 1] to uint8 and return a base64 string."""
    arr = np.clip(np.round(np.asarray(values, dtype=float) * 255.0), 0, 255).astype(np.uint8)
    return base64.b64encode(arr.tobytes()).decode("ascii")

def encode_pae_payload(pae_path):
    """Load a PAE matrix and pack it for the interactive report viewer.

    Returns None if the matrix is missing, empty or not square.
    """
    if is_missing_input(pae_path, "NO_FILE_PAE"):
        return None
    try:
        pae = np.genfromtxt(pae_path, delimiter="\t")
    except (ValueError, IOError):
        return None
    if pae.ndim != 2 or pae.shape[0] == 0 or pae.shape[0] != pae.shape[1]:
        return None
    if not np.isfinite(pae).all():
        pae = np.nan_to_num(pae, nan=PAE_MAX, posinf=PAE_MAX, neginf=0.0)
    scaled = np.clip(pae, 0.0, PAE_MAX) / PAE_MAX
    return {
        "size": int(pae.shape[0]),
        "max": PAE_MAX,
        "data": encode_u8_b64(scaled.ravel(order="C")),
    }

def get_residue_list(struct_file):
    """Ordered [chain id, residue number, residue name, insertion code, is_ligand]
    list for a structure.

    The order must match the residue/token order of the PAE matrix, which for
    prediction outputs follows the input sequence order: chains in file order,
    residues in ascending order within a chain. Waters are excluded because
    they are never scored tokens. HETATM residues (ligands, cofactors,
    modifications) are flagged so the report viewer can give them a fixed,
    selectable lane width on the PAE axes.
    """
    structure = parse_structure(struct_file, "report_ref")
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == "W":
                    continue
                hetflag, resseq, icode = residue.id
                residues.append([
                    chain.id,
                    int(resseq),
                    residue.resname,
                    icode.strip() or None,
                    hetflag != " ",
                ])
        break  # first model only - predictions carry a single model
    return residues

def load_msa_matrix(msa_path):
    """Parse an MSA and compute the sorted identity-scaled coverage matrix.

    Returns (n_cols, rows, identities) where each covered cell of a row
    holds that sequence's identity to the query (NaN for gaps), identities
    holds the per-row sequence identity to the query, and rows are ordered
    by ascending identity to the query. Returns None when the MSA file is
    missing.
    """
    if is_missing_input(msa_path):
        return None
    msa = []
    with open(msa_path, "r") as in_file:
        for line in in_file:
            msa.append([int(x) for x in line.strip().split()])
    if not msa:
        return None

    # Pad jagged MSAs to avoid shape errors in downstream plotting
    max_len = max(len(row) for row in msa)
    if any(len(row) != max_len for row in msa):
        msa = [row + [21] * (max_len - len(row)) for row in msa]

    seqid = []
    for sequence in msa:
        matches = [
            1.0 if first == other else 0.0 for first, other in zip(msa[0], sequence)
        ]
        seqid.append(sum(matches) / len(matches))

    seqid_sort = sorted(range(len(seqid)), key=seqid.__getitem__)

    non_gaps = []
    for sequence in msa:
        non_gaps.append(
            [float(num != 21) if num != 21 else float("nan") for num in sequence]
        )

    sorted_non_gaps = [non_gaps[i] for i in seqid_sort]
    final = []
    for sorted_seq, identity in zip(
        sorted_non_gaps, [seqid[i] for i in seqid_sort]
    ):
        final.append(
            [
                value * identity if not isinstance(value, str) else value
                for value in sorted_seq
            ]
        )

    return len(msa[0]), final, [seqid[i] for i in seqid_sort]

def build_msa_payload(msa_path):
    """Pack the sorted MSA coverage matrix for the interactive report viewer.

    Values are quantised to uint8. Deep MSAs are subsampled uniformly after the
    identity sort so the embedded report stays small; the full sequence count
    is reported alongside so the UI can flag the reduction.
    """
    loaded = load_msa_matrix(msa_path)
    if loaded is None:
        return None
    n_cols, rows, identities = loaded
    n_rows = len(rows)
    if n_rows > MSA_MAX_ROWS:
        keep = np.unique(np.round(np.linspace(0, n_rows - 1, MSA_MAX_ROWS)).astype(int))
        shown = [rows[i] for i in keep]
        shown_identities = [identities[i] for i in keep]
    else:
        shown = rows
        shown_identities = identities
    # Row 0 is the most similar sequence (the query): the original coverage
    # plot drew the identity-sorted rows with origin="lower", which puts the
    # query at the top of the image
    shown = shown[::-1]
    shown_identities = shown_identities[::-1]
    # NaN (gap) encodes to 0 exactly like a fully divergent covered column;
    # the per-row gap masks are implicit because gaps are rare and the raw
    # coverage overlay is derived from the matrix by the viewer.
    flat = np.array([
        0.0 if (isinstance(v, float) and np.isnan(v)) else v for row in shown for v in row
    ])
    seqid = shown_identities
    return {
        "nSeq": int(n_rows),
        "nSeqShown": int(len(shown)),
        "nCol": int(n_cols),
        "data": encode_u8_b64(flat),
        "seqid": encode_u8_b64(seqid),
    }

def is_missing_input(path, placeholder_prefix="NO_FILE"):
    return (
        not path
        or os.path.basename(path).startswith(placeholder_prefix)
        or not os.path.exists(path)
        or os.path.getsize(path) == 0
    )

def generate_msa_image(msa_path, name, out_dir, in_type):
    """Render the static MSA sequence-coverage PNG published with the run.

    The interactive report embeds the same data directly; the PNG stays as a
    lightweight artifact for downstream use.
    """
    loaded = load_msa_matrix(msa_path)
    if loaded is None:
        return
    xaxis_size, final, _ = loaded

    plt.figure(figsize=(16, 10), dpi=100)
    plt.title("Sequence coverage", fontsize=30, pad=36)
    plt.imshow(
        final,
        interpolation="nearest",
        aspect="auto",
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        origin="lower",
        extent=(0, xaxis_size, 0, len(final))
    )

    column_counts = [0] * xaxis_size
    for col in range(xaxis_size):
        for row in final:
            if not (isinstance(row[col], float) and np.isnan(row[col])):
                column_counts[col] += 1

    plt.plot(column_counts, color="black")
    plt.xlim(0, xaxis_size)
    plt.ylim(0, len(final))

    plt.tick_params(axis="both", which="both", labelsize=18)

    cbar = plt.colorbar()
    cbar.set_label("Sequence identity to query", fontsize=24, labelpad=24)
    cbar.ax.tick_params(labelsize=18)
    plt.xlabel("Positions", fontsize=24, labelpad=24)
    plt.ylabel("Sequences", fontsize=24, labelpad=36)
    plt.savefig(f"{out_dir}/{name}_{in_type}_seq_coverage.png")

def get_structure_parser(struct_file):
    suffix = os.path.splitext(struct_file)[1].lower()
    if suffix == ".pdb":
        return PDB.PDBParser(QUIET=True)
    if suffix in [".cif", ".mmcif"]:
        return PDB.MMCIFParser(QUIET=True)
    raise NotImplementedError("Reporting only supported for .pdb, .cif and .mmcif filetypes")


def parse_structure(struct_file, structure_id):
    parser = get_structure_parser(struct_file)
    return parser.get_structure(structure_id, struct_file)


def natural_path_sort_key(path):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(path))]


def pdb_to_lddt(struct_file, generate_tsv):
    """Per-residue pLDDT for the top-ranked structure, from B-factors.

    Returns (values, average). When generate_tsv is "y" the per-residue values
    are also written next to the structure file as <name>_plddt.tsv.
    """
    plddt_values = []

    structure = parse_structure(struct_file, "")

    for residue in structure.get_residues():
        res_pLDDT_tot = 0
        res_atom_count = 0

        for atom in residue.get_atoms():
            res_atom_count += 1
            res_pLDDT_tot += atom.get_bfactor()

        # Residue-level mean for ESMfold atom-level pLDDT
        res_pLDDT_ave = res_pLDDT_tot / res_atom_count

        if res_pLDDT_ave < 1.0:
            res_pLDDT_ave *= 100
        plddt_values.append(res_pLDDT_ave)

    average = round(sum(plddt_values) / len(plddt_values), 3) if plddt_values else 0.0

    if generate_tsv == "y":
        output_file = f"{os.path.splitext(struct_file)[0]}_plddt.tsv"
        with open(output_file, "w") as outfile:
            outfile.write(" ".join(map(str, plddt_values)) + "\n")

    return plddt_values, average


def read_ranked_score_tsv(tsv_path, model_count):
    scores = ["n/a"] * model_count
    if is_missing_input(tsv_path):
        return scores

    rows = []
    with open(tsv_path, "r") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            try:
                rank = int(row[0])
                score = float(row[1])
            except ValueError:
                continue
            rows.append((rank, f"{score:.3f}"))

    parsed_ranks = [rank for rank, _ in rows]
    rank_offset = 1 if parsed_ranks and 0 not in parsed_ranks and 1 in parsed_ranks else 0
    for rank, score in sorted(rows, key=lambda item: item[0]):
        rank -= rank_offset
        if 0 <= rank < model_count:
            scores[rank] = score

    return scores


def read_pair_score_tsv(tsv_path, model_count):
    scores = [{} for _ in range(model_count)]
    if is_missing_input(tsv_path):
        return scores

    with open(tsv_path, "r") as handle:
        reader = list(csv.reader(handle, delimiter="\t"))

    if not reader or len(reader[0]) < 2:
        return scores

    parsed_headers = []
    rank_headers = []
    for header in reader[0][1:]:
        try:
            parsed_rank = int(header)
            parsed_headers.append(parsed_rank)
            rank_headers.append(parsed_rank)
        except ValueError:
            rank_headers.append(None)

    rank_offset = 1 if parsed_headers and 0 not in parsed_headers and 1 in parsed_headers else 0

    for row in reader[1:]:
        if len(row) < 2:
            continue
        pair_label = row[0]
        for col_idx, value in enumerate(row[1:]):
            rank = rank_headers[col_idx] if col_idx < len(rank_headers) else None
            if rank is None:
                continue
            rank -= rank_offset
            if not (0 <= rank < model_count):
                continue
            try:
                parsed = float(value)
            except ValueError:
                continue
            # Native scores can be NaN for chains without an interface
            # (e.g. ligands); a NaN cell would render as an unreportable
            # block, so it is dropped here like any other missing value
            if not math.isfinite(parsed):
                continue
            scores[rank][pair_label] = f"{parsed:.4f}"

    return scores


def build_pair_score_matrices(score_maps):
    matrices = []

    for score_map in score_maps:
        if not score_map:
            matrices.append({"chains": [], "rows": []})
            continue

        chain_ids = sorted(
            {
                chain_id
                for pair_label in score_map.keys()
                for chain_id in pair_label.split(":")
                if ":" in pair_label
            }
        )
        rows = []
        for row_chain in chain_ids:
            row = []
            for col_chain in chain_ids:
                row.append(score_map.get(f"{row_chain}:{col_chain}", ""))
            rows.append(row)

        matrices.append({"chains": chain_ids, "rows": rows})

    return matrices


print("Starting...")

version = "1.0.0"
# Display names keyed by the shorthand mode tokens the pipeline passes as
# --type (the same tokens baked into the report file names); unknown tokens
# fall back to the token itself in programName (see below).
model_name_map = {
    "af2": "AlphaFold2",
    "af3": "AlphaFold3",
    "bz2": "Boltz2",
    "chai1": "Chai1",
    "esm2": "ESMFold2",
    "if2": "IntelliFold2",
    "ob0": "OpenBind0",
    "odde": "OpenDDE",
    "of3": "OpenFold3",
    "px2": "Protenix2",
    "rf3": "RosettaFold3",
}

parser = argparse.ArgumentParser()
parser.add_argument("--type", dest="in_type")
parser.add_argument("--tag", dest="file_tag", default=None, help="File-name tag; defaults to --type")
parser.add_argument(
    "--generate_tsv", choices=["y", "n"], default="n", dest="generate_tsv"
)
parser.add_argument("--msa", dest="msa", default=None)
parser.add_argument("--pdb", dest="pdb", required=True, nargs="+")
parser.add_argument("--pae", dest="pae", default=None)
parser.add_argument("--iptm", dest="iptm", default=None)
parser.add_argument("--ipsae", dest="ipsae", default=None)
parser.add_argument("--chainwise", dest="chainwise", default=None)
parser.add_argument("--chainwise_ipsae", dest="chainwise_ipsae", default=None)
parser.add_argument("--name", dest="name")
parser.add_argument("--output_dir", dest="output_dir")
parser.add_argument("--html_template", dest="html_template")
parser.add_argument("--version", action="version", version=f"{version}")
parser.set_defaults(output_dir="")
parser.set_defaults(in_type="af3")
parser.set_defaults(name="")
args = parser.parse_args()
file_tag = args.file_tag or args.in_type

# The report shows a single prediction per method: the top-ranked structure.
# Additional ranked models may be passed for provenance but are not embedded.
top_structure = sorted(args.pdb, key=natural_path_sort_key)[0]

_, plddt_average = pdb_to_lddt(top_structure, args.generate_tsv)

generate_msa_image(args.msa, args.name, args.output_dir, file_tag)

print("generating html report...")
iptm_score = read_ranked_score_tsv(args.iptm, 1)[0]
ipsae_score = read_ranked_score_tsv(args.ipsae, 1)[0]
# The chainwise iPTM store carries the per-chain pTM on the diagonal and
# the pairwise interface ipTM off it; the ipSAE store is a full
# asymmetric chain-pair matrix from colabfold's native interface scores.
chainwise_matrix = build_pair_score_matrices(
    read_pair_score_tsv(args.chainwise, 1)
)[0]
chainwise_ipsae_matrix = build_pair_score_matrices(
    read_pair_score_tsv(args.chainwise_ipsae, 1)
)[0]

# Residue labels follow the PAE/token order of the top-ranked model.
# Cross-highlighting between the PAE/MSA and the structure is only enabled
# when this list length matches the matrix dimensions (checked client-side).
residue_list = get_residue_list(top_structure)
pae_payload = encode_pae_payload(args.pae)
msa_payload = build_msa_payload(args.msa)

proteinfold_template = open(args.html_template, "r").read()

model_name = os.path.basename(top_structure)
model_format = "pdb" if os.path.splitext(model_name)[1].lower() == ".pdb" else "mmcif"
model_data = open(top_structure, "r").read()


def script_safe_json_dumps(obj):
    """json.dumps with <, > and & escaped, so a payload containing "</script>"
    cannot terminate the enclosing <script type="application/json"> element early."""
    return (
        json.dumps(obj)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


report_config = {
    "sampleName": args.name,
    "programName": model_name_map.get(args.in_type.lower(), args.in_type),
    "structFormat": model_format,
    "model": model_name,
    "modelData": model_data,
    "plddtAverage": plddt_average,
    "iptmScore": iptm_score,
    "ipsaeScore": ipsae_score,
    "chainwise": chainwise_matrix,
    "chainwise_ipsae": chainwise_ipsae_matrix,
    "residues": residue_list,
    "pae": pae_payload,
    "msa": msa_payload,
}
config_blob = (
    '<script type="application/json" id="report-config">'
    f"{script_safe_json_dumps(report_config)}</script>"
)
proteinfold_template = proteinfold_template.replace(
    "</head>", f"{config_blob}\n  </head>", 1
)

# The template marks optional sections with *-placeholder markers; the Python
# side swaps in mount points for data that exists and strips the whole section
# otherwise, so the page never shows empty panels.
if msa_payload is not None:
    proteinfold_template = proteinfold_template.replace(
        '<div id="seq_cov_placeholder"></div>',
        '<div id="msa-root" class="w-full"></div>',
    )
else:
    pattern = r'<section[^>]*id="seq_coverage_container".*?</section>'
    proteinfold_template = re.sub(pattern, "", proteinfold_template, flags=re.DOTALL)

if pae_payload is not None:
    proteinfold_template = proteinfold_template.replace(
        '<div id="pae_placeholder"></div>',
        '<div id="pae-root" class="w-full"></div>',
    )
else:
    pattern = r'<section[^>]*id="pae_container".*?</section>'
    proteinfold_template = re.sub(pattern, "", proteinfold_template, flags=re.DOTALL)


with open(
    f"{args.output_dir}/{args.name}_{file_tag.lower()}_report.html", "w"
) as out_file:
    out_file.write(proteinfold_template)
