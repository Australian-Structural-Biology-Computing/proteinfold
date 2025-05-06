from utils import reset_residue_numbers, align_structures, generate_plddt_plot, generate_sequence_coverage_plot
import os
import base64
from Bio import PDB
import argparse
from io import BytesIO

def generate_report(name, out_dir, structures, msa_files=None, type="standard", html_template=None):
    """
    Generate a report for protein structures, with optional comparison functionality.

    Args:
        name (str): Name of the report.
        out_dir (str): Output directory for the report.
        structures (list): List of structure file paths.
        msa_files (list, optional): List of MSA file paths. Defaults to None.
        comparison (bool, optional): Whether to generate a comparison report. Defaults to False.
        html_template (str, optional): Path to the HTML template for comparison. Defaults to None.
    """
    for structure in structures:
        if "esmfold" in structure:
            reset_residue_numbers(pdb_file, pdb_file) #Output pdb overwrite input to reset numbers

    if type == "comparison":
        aligned_structures = align_structures(structures)

        # Save the reference structure
        io = PDB.PDBIO()
        ref_structure_path = os.path.join(out_dir, "aligned_structure_0.pdb")
        io.set_structure(aligned_structures[0])
        io.save(ref_structure_path)
        aligned_structures[0] = ref_structure_path

        # Replace structures with aligned versions
        structures = aligned_structures

    template = open(html_template, "r").read()
    template = template.replace("*sample_name*", name)
    template = template.replace("*prog_name*", type)

    in_type = "NOT_COLABFOLD" # TODO: change args so that in_type can be distinguished from report type
    # Generate the sequence coverage plot and encode directly as base64 image rather than reading externally
    seq_cov_imgs = []
    if msa_files:
        for msa_file in msa_files:
            if msa_file and not msa_file.endswith("NO_FILE"):
                seq_cov_fig, _ = generate_sequence_coverage_plot(msa_file, out_dir, name, in_type, save_image=True)
                buf = BytesIO()
                seq_cov_fig.write_image(buf, format="png")
                buf.seek(0)
                encoded_image = base64.b64encode(buf.read()).decode("utf-8")
                seq_cov_img = f"data:image/png;base64,{encoded_image}"
                seq_cov_imgs.append(seq_cov_img)

    if type == "comparison":
        args_msa_array_js = (f"""const SEQ_COV_IMGS = [{", ".join([f'"{fig}"' for fig in seq_cov_imgs])}];""")
        template = template.replace("const SEQ_COV_IMGS = [];", args_msa_array_js)
    else:
        template = template.replace("seq_coverage.png", f"data:image/png;base64,{seq_cov_imgs[0]}")

    # Generate the pLDDT plot
    fig = generate_plddt_plot(structures)
    html_content = fig.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
    )
    with open(f"{out_dir}/{name}_coverage_LDDT.html", "w") as out_file:
        out_file.write(html_content)

    if type == "comparison":
        # Add aligned structures
        args_pdb_array_js = (
            "const MODELS = [" + ",\n".join([f'"{model}"' for model in structures]) + "];"
        )
        template = template.replace("const MODELS = [];", args_pdb_array_js)

    with open(f"{out_dir}/{name}_coverage_LDDT.html", "r") as lddt_file:
        lddt_html = lddt_file.read()
        template = template.replace('<div id="lddt_placeholder"></div>', lddt_html)

    with open(f"{out_dir}/{name}_{type}_report.html", "w") as out_file:
        out_file.write(template)

def main():
    parser = argparse.ArgumentParser(description="Generate protein structure reports.")
    parser.add_argument("--name", required=True, help="Name of the report.")
    parser.add_argument("--output_dir", required=True, help="Output directory for the report.")
    parser.add_argument("--structs", required=True, nargs="+", help="List of structure file paths.")
    parser.add_argument("--msa", nargs="+", default=None, help="List of MSA file paths (optional).")
    parser.add_argument("--type", required=False, default="standard", choices=["standard", "comparison"], help="The type of report file generated .") # TODO: change to --type with options in case there are other reports
    parser.add_argument("--html_template", default=None, help="Path to the HTML template for comparison (optional).")

    args = parser.parse_args()

    print("Generating report.....")

    generate_report(
        name=args.name,
        out_dir=args.output_dir,
        structures=args.structs,
        msa_files=args.msa,
        type=args.type,
        html_template=args.html_template,
    )

if __name__ == "__main__":
    main()
