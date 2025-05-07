from utils import reset_residue_numbers, align_structures, plddt_from_struct_b_factor, generate_plddt_plot, generate_pae_plot, generate_sequence_coverage_plot
import os
import base64
from Bio import PDB
import argparse
from io import BytesIO

def generate_report(name, out_dir, structures, msa_files=None, pae_files=None, type="standard", html_template=None, write_htmls=True, seq_cov_as_html=False):
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

    print("Structures:", structures)

    template = open(html_template, "r").read()
    template = template.replace("*sample_name*", name)
    # TODO: pass in program name as an argument
    #template = template.replace("*prog_name*", prog_name)

    #Implement an LDDT averages function
    #plddt_rank_0 = plddt_from_struct_b_factor(structures[0])
    #template = template.replace("*prog_name*", prog_name)

    # Populate MODELS into the HTML template
    models_js = ("const MODELS = [" + ",\n".join([f'"{model}"' for model in structures]) + "];")
    template = template.replace("const MODELS = [];", models_js)
    # Populate MODELS_DATA with the content of the PDB files
    pdb_strings = [open(structure, "r").read().replace("\n", "\\n") for structure in structures]
    models_data = ",\n".join([f'"{pdb_string}"' for pdb_string in pdb_strings])
    models_data_js = f"const MODELS_DATA = [{models_data}];"
    template = template.replace("const MODELS_DATA = [];", models_data_js)

    # Generate sequence coverage plots and convert to HTML
    if msa_files:
        for msa_file in msa_files:
            if msa_file and not msa_file.endswith("NO_FILE"):
                seq_cov_fig, seq_cov_img_path = generate_sequence_coverage_plot(msa_file, out_dir, name, in_type="NOT_COLABFOLD", save_image=True)
                seq_cov_img_encoded = base64.b64encode(open(seq_cov_img_path, "rb").read()).decode("utf-8")
                seq_cov_img_tag = f'<img src="data:image/png;base64,{seq_cov_img_encoded}" alt="Sequence Coverage Image">'

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

    if pae_files:
        pae_figs = []
        for pae_file in pae_files:
            pae_figs.append(generate_pae_plot(pae_file, out_dir, name, save_image=True))
        pae_html = pae_figs[0].to_html(
            full_html=False,
            include_plotlyjs="cdn",
            config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
        )
   #Generate PAE plot and conver to HTML TODO: currently onlt the first

    # Place the HTML plots in their div conttainer
    if seq_cov_as_html == True:
        template = template.replace('<div id="seq_cov_placeholder"></div>', seq_cov_html)
    else:
        template = template.replace('<div id="seq_cov_placeholder"></div>', seq_cov_img_tag)
    template = template.replace('<div id="lddt_placeholder"></div>', plddt_html)
    template = template.replace('<div id="pae_placeholder"></div>', pae_html)

    if write_htmls:
        with open(f"{out_dir}/{name}_coverage_pLDDT.html", "w") as out_file:
            out_file.write(seq_cov_html)
        with open(f"{out_dir}/{name}_coverage_MSA.html", "w") as out_file:
            out_file.write(seq_cov_html)

    # if type == "comparison":
      #  args_msa_array_js = (f"""const SEQ_COV_IMGS = [{", ".join([f'"{fig}"' for fig in seq_cov_figs])}];""")
      #  template = template.replace("const SEQ_COV_IMGS = [];", args_msa_array_js)

    with open(f"{out_dir}/{name}_{type}_report.html", "w") as out_file:
        out_file.write(template)

def main():
    parser = argparse.ArgumentParser(description="Generate protein structure reports.")
    parser.add_argument("--name", required=True, help="Name of the report.")
    parser.add_argument("--output_dir", required=True, help="Output directory for the report.")
    parser.add_argument("--structs", required=True, nargs="+", help="List of structure file paths.")
    parser.add_argument("--msa", nargs="+", default=None, help="List of MSA file paths (optional).")
    parser.add_argument("--pae", nargs="+", default=None, help="List of PAE file paths (optional).")
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
        pae_files=args.pae,
        type=args.type,
        html_template=args.html_template,
        write_htmls=args.write_htmls,
        seq_cov_as_html=False,
    )

if __name__ == "__main__":
    main()
