#!/usr/bin/env python

import math
import os
import sys
import argparse
import glob
import json
from matplotlib import pyplot as plt
import numpy as np
from collections import OrderedDict
import base64
import plotly.graph_objects as go
import re
from Bio import PDB

def is_missing_input(path):
    """True when a metric input is absent or empty.

    The inference task signals an unavailable MSA (the ESMFold2
    language-model variants fold without one) with a zero-size TSV, which
    this treats the same as a missing file.
    """
    return not path or not os.path.exists(path) or os.path.getsize(path) == 0

def generate_pae_plot(pae_matrix, out_dir, name, save_image=False):
    #save_image=False because plotly needs a local install of Google Chrome to save images.....
    """
    Generate a Plotly heatmap for Predicted Aligned Error (PAE) data.

    Args:
        pae_matrix (2D array-like): The PAE matrix (from the packed metric store).
    Returns:
        fig: A Plotly figure object of the PAE heatmap in green color scale
    """
    pae = np.asarray(pae_matrix, dtype=float)
    fig = go.Figure()

    # Add heatmap
    fig.add_trace(
        go.Heatmap(
            z=pae,
            colorscale="Greens_r",
            zmin=0,
            zmax=30,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Scored Residue", minallowed=0, maxallowed=pae.shape[0]-1),
        yaxis=dict(title="Aligned Residue", minallowed=0, maxallowed=pae.shape[1]-1, autorange="reversed"),
        width=600,
        height=600,
    )

    if save_image:
            image_path = f"{out_dir}/{name+('_' if name else '')}PAE.png"
            fig.write_image(image_path, width=800, height=800)

    return fig

def generate_output_images(msa_path, plddt_data, name, out_dir, in_type, generate_tsv, pdb, pae_matrix=None):
    msa = []
    if not is_missing_input(msa_path):
        if str(msa_path).endswith(".a3m"):
            # a3m MSA (packed-metrics mode): integer-encode like extract_metrics
            sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
            from extract_metrics import a3m_to_int

            msa = [list(row) for row in a3m_to_int(msa_path)]
        else:
            with open(msa_path, "r") as in_file:
                for line in in_file:
                    msa.append([int(x) for x in line.strip().split()])

        # Pad jagged MSAs to avoid shape errors in downstream plotting
        if msa:
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

        xaxis_size = len(final[0])
        yaxis_size = len(final)

        # ##################################################################
        plt.figure(figsize=(16, 10), dpi=100)
        # ##################################################################
        plt.title("Sequence coverage", fontsize=30, pad=36)
        plt.imshow(
            final,
            interpolation="nearest",
            aspect="auto",
            cmap="rainbow_r",
            vmin=0,
            vmax=1,
            origin="lower",
            extent=(0, xaxis_size, 0, yaxis_size)
        )

        column_counts = [0] * len(msa[0])
        for col in range(len(msa[0])):
            for row in msa:
                if row[col] != 21:
                    column_counts[col] += 1

        plt.plot(column_counts, color="black")
        plt.xlim(0, len(msa[0]))
        plt.ylim(0, len(msa))

        plt.tick_params(axis="both", which="both", labelsize=18)

        cbar = plt.colorbar()
        cbar.set_label("Sequence identity to query", fontsize=24, labelpad=24)
        cbar.ax.tick_params(labelsize=18)
        plt.xlabel("Positions", fontsize=24, labelpad=24)
        plt.ylabel("Sequences", fontsize=24, labelpad=36)
        plt.savefig(f"{out_dir}/{name}_{in_type}_seq_coverage.png")

        # ##################################################################

    plddt_per_model = OrderedDict()
    output_data = plddt_data

    if generate_tsv == "y":
        for plddt_path in output_data:
            with open(plddt_path, "r") as in_file:
                plddt_per_model[os.path.basename(plddt_path)[:-4]] = [
                    float(x) for x in in_file.read().strip().split()
                ]
    else:
        for i, plddt_values_str in enumerate(output_data):
            plddt_per_model[i] = []
            plddt_per_model[i] = [float(x) for x in plddt_values_str.strip().split()]

    fig = go.Figure()
    for idx, (model_name, value_plddt) in enumerate(plddt_per_model.items()):
        rank_label = os.path.splitext(pdb[idx])[0]
        fig.add_trace(
            go.Scatter(
                x=list(range(len(value_plddt))),
                y=value_plddt,
                mode="lines",
                name=rank_label,
                text=[f"({i}, {value:.2f})" for i, value in enumerate(value_plddt)],
                hoverinfo="text",
            )
        )
    fig.update_layout(
        title=dict(text="Predicted LDDT per position", x=0.5, xanchor="center"),
        xaxis=dict(
            title="Positions", showline=True, linecolor="black", gridcolor="WhiteSmoke", minallowed=0, maxallowed=len(value_plddt)-1
        ),
        yaxis=dict(
            title="Predicted LDDT",
            range=[0, 100],
            fixedrange=True,
            showline=True,
            linecolor="black",
            gridcolor="WhiteSmoke",
        ),
        legend=dict(yanchor="bottom", y=0.02, xanchor="right", x=1, bordercolor="Black", borderwidth=1),
        plot_bgcolor="white",
        width=600,
        height=600,
        modebar_remove=["toImage", "zoomIn", "zoomOut"],
    )
    html_content = fig.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
    )

    with open(
        f"{out_dir}/{name+('_' if name else '')}coverage_LDDT.html", "w"
    ) as out_file:
        out_file.write(html_content)

    if pae_matrix is not None:
        pae_fig = generate_pae_plot(pae_matrix, out_dir, name)
        pae_html_content = pae_fig.to_html(
            full_html=False,
            include_plotlyjs="cdn",
            config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
        )
        with open(
            f"{out_dir}/{name+('_' if name else '')}PAE.html", "w"
        ) as pae_out_file:
            pae_out_file.write(pae_html_content)

def generate_plots(msa_path, plddt_paths, name, out_dir):
    msa = []
    with open(msa_path, "r") as in_file:
        for line in in_file:
            msa.append([int(x) for x in line.strip().split()])

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
    for sorted_seq, identity in zip(sorted_non_gaps, [seqid[i] for i in seqid_sort]):
        final.append(
            [
                value * identity if not isinstance(value, str) else value
                for value in sorted_seq
            ]
        )

    # Plotting Sequence Coverage using Plotly
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=final,
            colorscale="Rainbow",
            zmin=0,
            zmax=1,
        )
    )
    fig.update_layout(
        title="Sequence coverage", xaxis_title="Positions", yaxis_title="Sequences"
    )
    # Save as interactive HTML instead of an image
    fig.savefig(f"{out_dir}/{name+('_' if name else '')}seq_coverage.png")

    # Plotting Predicted LDDT per position using Plotly
    plddt_per_model = OrderedDict()
    plddt_paths.sort()
    for plddt_path in plddt_paths:
        with open(plddt_path, "r") as in_file:
            plddt_per_model[os.path.basename(plddt_path)[:-4]] = [
                float(x) for x in in_file.read().strip().split()
            ]

    i = 0
    for model_name, value_plddt in plddt_per_model.items():
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=list(range(len(value_plddt))),
                y=value_plddt,
                mode="lines",
                name=model_name,
            )
        )
        fig.update_layout(title="Predicted LDDT per Position")
        fig.savefig(f"{out_dir}/{name+('_' if name else '')}coverage_LDDT_{i}.png")
        i += 1



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


def get_atom_id(atom):
    residue = atom.get_parent()
    chain = residue.get_parent()
    return (chain.get_id(), residue.get_id(), atom.name)


def get_non_hydrogen_atoms_by_id(structure):
    return {
        get_atom_id(atom): atom
        for atom in structure.get_atoms()
        if atom.element != "H"
    }


def align_structures(structures):
    parsed_structures = [
        parse_structure(structure, f"Structure_{i}")
        for i, structure in enumerate(structures)
    ]
    ref_structure = parsed_structures[0]
    atom_maps = [get_non_hydrogen_atoms_by_id(structure) for structure in parsed_structures]
    common_atom_ids = set(atom_maps[0])

    for atom_map in atom_maps[1:]:
        common_atom_ids.intersection_update(atom_map)

    if not common_atom_ids:
        raise ValueError("No common non-hydrogen atoms found between structures.")

    ref_atom_ids = [atom_id for atom_id in atom_maps[0] if atom_id in common_atom_ids]
    ref_atoms = [atom_maps[0][atom_id] for atom_id in ref_atom_ids]
    super_imposer = PDB.Superimposer()
    aligned_structures = [ref_structure]

    for i, structure in enumerate(parsed_structures[1:], start=1):
        target_atoms = [atom_maps[i][atom_id] for atom_id in ref_atom_ids]
        super_imposer.set_atoms(ref_atoms, target_atoms)
        super_imposer.apply(structure.get_atoms())

        aligned_structure = f"aligned_structure_{i}.cif"
        io = PDB.MMCIFIO()
        io.set_structure(structure)
        io.save(aligned_structure)
        aligned_structures.append(aligned_structure)

    return aligned_structures


def natural_path_sort_key(path):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(path))]


def pdb_to_lddt(struct_files, generate_tsv):
    struct_files_sorted = sorted(struct_files, key=natural_path_sort_key)

    output_lddt = []
    averages = []

    for struct_file in struct_files_sorted:
        plddt_values = []

        structure = parse_structure(struct_file, "")

        for residue in structure.get_residues():
            res_pLDDT_tot = 0
            res_atom_count = 0

            for atom in residue.get_atoms():
                res_atom_count +=1
                res_pLDDT_tot += atom.get_bfactor()

            # Residue-level mean for ESMfold atom-level pLDDT
            res_pLDDT_ave = res_pLDDT_tot/res_atom_count

            if res_pLDDT_ave < 1.0:
                res_pLDDT_ave *= 100
            plddt_values.append(res_pLDDT_ave)

        # Calculate the average PLDDT value for the current file
        if plddt_values:
            avg_plddt = sum(plddt_values) / len(plddt_values)
            averages.append(round(avg_plddt, 3))
        else:
            averages.append(0.0)

        if generate_tsv == "y":
            output_file = f"{os.path.splitext(struct_file)[0]}_plddt.tsv"
            with open(output_file, "w") as outfile:
                outfile.write(" ".join(map(str, plddt_values)) + "\n")
            output_lddt.append(output_file)
        else:
            plddt_values_string = " ".join(map(str, plddt_values))
            output_lddt.append(plddt_values_string)

    return output_lddt, averages


def read_ranked_scores_from_store(scalar_records, sample_id, metric, model_count):
    """Per-rank scores for one sample, sourced from the packed scalar store.

    Semantics identical to read_ranked_score_tsv: values are rendered at 3
    decimals and ranks are shifted to 0-based positions when the native
    ranking starts at 1.
    """
    scores = ["n/a"] * model_count
    per_rank = {}
    for record in scalar_records:
        if record["sample_id"] != sample_id or record["metric"] != metric or record.get("chain_pair"):
            continue
        try:
            per_rank[int(record["rank"])] = f"{float(record['value']):.3f}"
        except (TypeError, ValueError):
            continue
    rank_offset = 1 if per_rank and 0 not in per_rank and 1 in per_rank else 0
    for rank in sorted(per_rank):
        index = rank - rank_offset
        if 0 <= index < model_count:
            scores[index] = per_rank[rank]
    return scores


def read_pair_scores_from_store(scalar_records, sample_id, metric, model_count):
    """Per-rank chain-pair score maps for one sample, from the packed store.

    Semantics identical to read_pair_score_tsv.
    """
    scores = [{} for _ in range(model_count)]
    per_rank = {}
    for record in scalar_records:
        if record["sample_id"] != sample_id or record["metric"] != metric or not record.get("chain_pair"):
            continue
        try:
            parsed = float(record["value"])
        except (TypeError, ValueError):
            continue
        # Native scores can be NaN for chains without an interface
        # (e.g. ligands); a NaN cell would render as an unreportable
        # block, so it is dropped here like any other missing value
        if not math.isfinite(parsed):
            continue
        per_rank.setdefault(int(record["rank"]), {})[record["chain_pair"]] = f"{parsed:.4f}"
    rank_offset = 1 if per_rank and 0 not in per_rank and 1 in per_rank else 0
    for rank, pairs in per_rank.items():
        index = rank - rank_offset
        if 0 <= index < model_count:
            scores[index] = pairs
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
model_name = {
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

def render_sample_report(name, in_type, struct_files, msa_path, scalar_records, array_store,
                         html_template, output_dir, generate_tsv="n", file_tag=None):
    """Render the report artefacts for one sample from packed metric stores."""
    file_tag = file_tag or in_type
    structures = sorted(struct_files, key=natural_path_sort_key)
    # Legacy reports plotted the rank-1 PAE (the only per-sample PAE the
    # pipeline emitted); keep that behaviour when reading from the store.
    pae_matrix = array_store.get(name, 1, "pae")
    lddt_data, lddt_averages = pdb_to_lddt(struct_files, generate_tsv)

    generate_output_images(
        msa_path, lddt_data, name, output_dir, in_type, generate_tsv, struct_files,
        pae_matrix=pae_matrix
    )

    print("generating html report...")
    iptm_scores = read_ranked_scores_from_store(scalar_records, name, "iptm", len(structures))
    ipsae_scores = read_ranked_scores_from_store(scalar_records, name, "ipsae", len(structures))
    chainwise_iptm_scores = read_pair_scores_from_store(scalar_records, name, "chainwise_iptm", len(structures))
    chainwise_ipsae_scores = read_pair_scores_from_store(scalar_records, name, "chainwise_ipsae", len(structures))
    chainwise_iptm_matrices = build_pair_score_matrices(chainwise_iptm_scores)
    chainwise_ipsae_matrices = build_pair_score_matrices(chainwise_ipsae_scores)
    aligned_structures = align_structures(structures)

    io = PDB.MMCIFIO()
    ref_structure_path = "aligned_structure_0.cif"
    io.set_structure(aligned_structures[0])
    io.save(ref_structure_path)
    aligned_structures[0] = ref_structure_path

    proteinfold_template = open(html_template, "r").read()

    model_names = [
        f"{os.path.splitext(os.path.basename(model))[0]}.cif" for model in structures
    ]
    models_data = [open(s, "r").read() for s in aligned_structures]


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
        "reportType": "standard",
        "sampleName": name,
        "programName": model_name.get(in_type.lower(), in_type),
        "structFormat": "cif",
        "models": model_names,
        "models_data": models_data,
        "lddt_averages": lddt_averages,
        "iptm_scores": iptm_scores,
        "ipsae_scores": ipsae_scores,
        "chainwise": chainwise_iptm_matrices,
        "chainwise_ipsae": chainwise_ipsae_matrices,
    }
    config_blob = (
        '<script type="application/json" id="report-config">'
        f"{script_safe_json_dumps(report_config)}</script>"
    )
    proteinfold_template = proteinfold_template.replace(
        "</head>", f"{config_blob}\n  </head>", 1
    )

    if not is_missing_input(msa_path):
        image_path = f"{output_dir}/{name}_{file_tag}_seq_coverage.png"
        with open(image_path, "rb") as in_file:
            data_uri = f"data:image/png;base64,{base64.b64encode(in_file.read()).decode('utf-8')}"
            proteinfold_template = proteinfold_template.replace(
                '<div id="seq_cov_placeholder"></div>',
                f'<img src="{data_uri}" alt="Sequence coverage (MSA)" '
                'class="w-full h-auto rounded" />',
            )
    else:
        pattern = r'<div id="seq_coverage_container".*?>.*?(<!--.*?-->.*?)*?</div>\s*</div>\s*</div>\s*</div>'
        proteinfold_template = re.sub(pattern, "", proteinfold_template, flags=re.DOTALL)

    with open(
        f"{output_dir}/{name + ('_' if name else '')}coverage_LDDT.html",
        "r",
    ) as in_file:
        lddt_html = in_file.read()
        proteinfold_template = proteinfold_template.replace(
            '<div id="lddt_placeholder"></div>', lddt_html
        )

    if pae_matrix is not None:
        with open(
            f"{output_dir}/{name + ('_' if name else '')}PAE.html",
            "r",
        ) as pae_in_file:
            pae_html = pae_in_file.read()
            proteinfold_template = proteinfold_template.replace(
                '<div id="pae_placeholder"></div>', pae_html
            )
    else:
        pattern = r'<div id="pae_container".*?>.*?(<!--.*?-->.*?)*?</div>\s*</div>'
        proteinfold_template = re.sub(pattern, "", proteinfold_template, flags=re.DOTALL)

    with open(
        f"{output_dir}/{name}_{file_tag.lower()}_report.html", "w"
    ) as out_file:
        out_file.write(proteinfold_template)

    # The coverage/PAE visualisations above are embedded into the report
    # itself; remove the transient files so they do not accumulate in the
    # task directory (only the report HTML is published).
    for transient in (
        f"{output_dir}/{name}_{file_tag}_seq_coverage.png",
        f"{output_dir}/{name}{'_' if name else ''}coverage_LDDT.html",
        f"{output_dir}/{name}_PAE.html",
        "aligned_structure_0.cif",
    ):
        try:
            os.remove(transient)
        except FileNotFoundError:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", dest="in_type")
    parser.add_argument(
        "--generate_tsv", choices=["y", "n"], default="n", dest="generate_tsv"
    )
    parser.add_argument("--manifest", required=True,
                        help="TSV: sample_id, packed scalars store, packed arrays store")
    parser.add_argument("--struct-ext", dest="struct_ext", default="cif", choices=["cif", "pdb"])
    parser.add_argument("--output_dir", dest="output_dir")
    parser.add_argument("--html_template", dest="html_template")
    parser.add_argument("--version", action="version", version=f"{version}")
    parser.set_defaults(output_dir="")
    parser.set_defaults(in_type="af3")
    args = parser.parse_args()

    sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
    import extract_metrics
    import metrics_store

    with open(args.manifest) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    for fields in rows:
        if len(fields) != 4:
            sys.exit(f"generate_report: malformed manifest line: {fields}")
        sample_id, scalars_path, arrays_path, file_tag = fields

        struct_files = extract_metrics.sort_paths_by_rank(
            glob.glob(f"{sample_id}_rank_*.{args.struct_ext}")
        )
        if not struct_files:
            sys.exit(f"generate_report: no ranked structures staged for {sample_id}")

        scalar_records = list(metrics_store.read_scalars(scalars_path))
        with metrics_store.ArrayStore(arrays_path) as array_store:
            render_sample_report(
                name=sample_id,
                in_type=args.in_type,
                file_tag=file_tag,
                struct_files=struct_files,
                msa_path=f"{sample_id}.a3m",
                scalar_records=scalar_records,
                array_store=array_store,
                html_template=args.html_template,
                output_dir=args.output_dir,
                generate_tsv=args.generate_tsv,
            )
        print(f"generate_report: rendered report for {sample_id}")

    print(f"generate_report: rendered {len(rows)} report(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
