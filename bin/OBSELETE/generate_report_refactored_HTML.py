from utils import reset_residue_numbers, align_structures, generate_plddt_plot, generate_sequence_coverage_plot
import os
import base64
from Bio import PDB
import argparse
from io import BytesIO

def generate_report(name, out_dir, structures, msa_files=None, type="standard", html_template=None, write_htmls=True):
    for structure in structures:
        if "esmfold" in structure:
            reset_residue_numbers(pdb_file, pdb_file) #Output pdb overwrite input to reset numbers

    template = open(html_template, "r").read()
    template = template.replace("*sample_name*", name)
    template = template.replace("*prog_name*", type)

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

    # Add aligned structures to report
    if type == "comparison":
        args_pdb_array_js = (
            "const MODELS = [" + ",\n".join([f'"{model}"' for model in structures]) + "];"
        )
        template = template.replace("const MODELS = [];", args_pdb_array_js)

    # Generate sequence coverage plots and convert to HTML
    if msa_files:
        for msa_file in msa_files:
            if msa_file and not msa_file.endswith("NO_FILE"):
                seq_cov_fig = generate_sequence_coverage_plot(msa_file, out_dir, name, in_type="NOT_COLABFOLD", save_image=True)
                seq_cov_html = seq_cov_fig.to_html(
                    full_html=False,
                    include_plotlyjs="cdn",
                    config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
                )

    # Generate the pLDDT plot and convert to HTML
    plddt_fig = generate_plddt_plot(structures)
    plddt_html = plddt_fig.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
    )

    # Place the HTML plots in their fdiv conttainer
    template = template.replace('<div id="seq_cov_placeholder"></div>', seq_cov_html)
    template = template.replace('<div id="lddt_placeholder"></div>', plddt_html)
    if write_htmls:
        with open(f"{out_dir}/{name}_coverage_pLDDT.html", "w") as out_file:
            out_file.write(seq_cov_html)
        with open(f"{out_dir}/{name}_coverage_MSA.html", "w") as out_file:
            out_file.write(seq_cov_html)

    if type == "comparison":
        args_msa_array_js = (f"""const SEQ_COV_IMGS = [{", ".join([f'"{fig}"' for fig in seq_cov_figs])}];""")
        template = template.replace("const SEQ_COV_IMGS = [];", args_msa_array_js)

    with open(f"{out_dir}/{name}_{type}_report.html", "w") as out_file:
        out_file.write(template)

def main():
    parser = argparse.ArgumentParser(description="Generate protein structure reports.")
    parser.add_argument("--name", required=True, help="Name of the report.")
    parser.add_argument("--output_dir", required=True, help="Output directory for the report.")
    parser.add_argument("--structs", required=True, nargs="+", help="List of structure file paths.")
    parser.add_argument("--msa", nargs="+", default=None, help="List of MSA file paths (optional).")
    parser.add_argument("--type", default="standard", choices=["standard", "comparison"], help="The type of report file generated .") # TODO: change to --type with options in case there are other reports
    parser.add_argument("--html_template", default=None, help="Path to the HTML template for comparison (optional).")
    parser.add_argument("--write_htmls", default=True, help="Write out seperate files for each html plot (optional).")

    args = parser.parse_args()

    print("Generating report.....")

    generate_report(
        name=args.name,
        out_dir=args.output_dir,
        structures=args.structs,
        msa_files=args.msa,
        type=args.type,
        html_template=args.html_template,
        write_htmls=args.write_htmls,
    )

if __name__ == "__main__":
    main()
