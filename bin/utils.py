from collections import OrderedDict
import plotly.graph_objects as go
from Bio import PDB
import matplotlib.pyplot as plt
import numpy as np

def reset_residue_numbers(input_pdb, output_pdb):  #TODO: use PDBIO instead of file I/O
    """
    Resets residue numbers (column 23-26) in a PDB file so the position starts from 1 for each chain
    and increment only when encountering a new residue.
    """
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        current_residue_number = 1
        previous_residue_id = None
        previous_chain = None

        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]  # Extract the chain identifier (column 22)
                residue_id = line[22:26].strip()  # Extract the residue ID (column 23-26)

                # Reset residue numbering if the chain changes
                if chain != previous_chain:
                    current_residue_number = 1
                    previous_chain = chain
                    previous_residue_id = None

                # Increment residue number if it's a new residue
                if residue_id != previous_residue_id:
                    if previous_residue_id is not None:  # Only increment after the first residue
                        current_residue_number += 1
                    previous_residue_id = residue_id

                # Update the line with the new residue number
                updated_line = (
                    line[:22] +
                    f"{current_residue_number:4}" +
                    line[26:]
                )
                outfile.write(updated_line)

            else:
                # Write non-ATOM/HETATM lines (e.g., TER, PARENT) without changes
                outfile.write(line)

def align_structures(structures):
    parser = PDB.PDBParser(QUIET=True)  #TODO: or use MMCIFParser for mmCIF files
    structures = [
        parser.get_structure(f"Structure_{i}", pdb) for i, pdb in enumerate(structures)
    ]
    ref_structure = structures[0]

    common_atoms = set(
        f"{atom.get_parent().get_id()[1]}-{atom.name}"
        for atom in ref_structure.get_atoms()
    )
    for i, structure in enumerate(structures[1:], start=1):
        common_atoms = common_atoms.intersection(
            set(
                f"{atom.get_parent().get_id()[1]}-{atom.name}"
                for atom in structure.get_atoms()
            )
        )

    ref_atoms = [
        atom
        for atom in ref_structure.get_atoms()
        if f"{atom.get_parent().get_id()[1]}-{atom.name}" in common_atoms
    ]
    # print(ref_atoms)
    super_imposer = PDB.Superimposer()
    aligned_structures = [structures[0]]  # Include the reference structure in the list

    for i, structure in enumerate(structures[1:], start=1):
        target_atoms = [
            atom
            for atom in structure.get_atoms()
            if f"{atom.get_parent().get_id()[1]}-{atom.name}" in common_atoms
        ]

        super_imposer.set_atoms(ref_atoms, target_atoms)
        super_imposer.apply(structure.get_atoms())

        aligned_structure = f"aligned_structure_{i}.pdb"
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(aligned_structure)
        aligned_structures.append(aligned_structure)

    return aligned_structures

def plddt_from_struct_b_factor(struct_file):
    """
    Uses the BioPython PDB package to extract residue plddt values from the b-factor column. Iterates over PDB objects rather than processes raw file
    """
    if str(struct_file).endswith(".pdb"):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(id=id, file=struct_file)
    elif str(struct_file).endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure(structure_id=id, filename=struct_file)
    else:
        print(f"{struct_file} is neither a PDB or mmCIF file!")

    res_list = []
    res_plddts = []
    plddt_tot = 0

    for model in structure:
        for chain in model:
            chain_res_list = chain.get_unpacked_list()
            res_list.extend(chain_res_list)
            for residue in chain:
                atom_list = residue.get_unpacked_list()
                atom_plddt_tot = 0
                for atom in residue:  # ESMFold and others have separate atom-wise values, so doing atom-wise to cover that and residue-wise
                    atom_plddt = atom.get_bfactor()
                    atom_plddt_tot += atom_plddt

                res_plddt = float(atom_plddt_tot / len(atom_list))

                if (res_plddt < 1):  # RFAA the multiplication of mean isn't failing. Anyway covering to a [0,100] range for any structure file1
                    res_plddt *= 100

                res_plddts.append(res_plddt)
                plddt_tot += res_plddt

    res_plddts = np.array(res_plddts)
    res_plddts = np.round(res_plddts, 2)

    return res_plddts

def generate_plddt_plot(structures):
    """
    Generate a Plotly figure for predicted lddt per position for given structures.

    Args:
        structures (list): List of structure file paths.

    Returns:
        go.Figure: Plotly figure object with lddt data.
    """
    plddt_per_struct = OrderedDict()

    # Extract plddt values for each structure
    for struct in structures:
        plddt_per_struct[struct] = plddt_from_struct_b_factor(struct)

    # Create the Plotly figure
    fig = go.Figure()

    for idx, (struct, plddts) in enumerate(plddt_per_struct.items()):
        fig.add_trace(
            go.Scatter(
                x=list(range(len(plddts))),
                y=plddts,
                mode="lines",
                name=f"rank-{idx}",
                text=[f"({idx}, {value:.2f})" for idx, value in enumerate(plddts)],
                hoverinfo="text",
            )
        )

    # Update layout
    fig.update_layout(
        title=dict(text="Predicted lddt per position", x=0.5, xanchor="center"),
        xaxis=dict(
            title="Positions", showline=True, linecolor="black", gridcolor="WhiteSmoke"
        ),
        yaxis=dict(
            title="Predicted lddt",
            range=[0, 100],
            showline=True,
            linecolor="black",
            gridcolor="WhiteSmoke",
        ),
        legend=dict(
            yanchor="bottom", y=0.02, xanchor="right", x=1, bordercolor="Black", borderwidth=1
        ),
        plot_bgcolor="white",
        width=600,
        height=600,
    )

    return fig

def generate_sequence_coverage_plot(msa_path, out_dir, name, in_type="standard", save_image=True):
    msa = []
    if in_type.lower() != "colabfold" and not msa_path.endswith("NO_FILE"):
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
        for sorted_seq, identity in zip(
            sorted_non_gaps, [seqid[i] for i in seqid_sort]
        ):
            final.append(
                [
                    value * identity if not isinstance(value, str) else value
                    for value in sorted_seq
                ]
            )

        # TODO: don't have a seperate save iamge plot and a plotly ploy
        # Plot the sequence coverage and save as image
        # ##################################################################
        if save_image:
            image_path = f"{out_dir}/{name+('_' if name else '')}seq_coverage.png"
            plt.figure(figsize=(14, 14), dpi=100)
            plt.title("Sequence coverage", fontsize=30, pad=36)
            plt.imshow(
                final,
                interpolation="nearest",
                aspect="auto",
                cmap="rainbow_r",
                vmin=0,
                vmax=1,
                origin="lower",
            )

            column_counts = [0] * len(msa[0])
            for col in range(len(msa[0])):
                for row in msa:
                    if row[col] != 21:
                        column_counts[col] += 1

            plt.plot(column_counts, color="black")
            plt.xlim(-0.5, len(msa[0]) - 0.5)
            plt.ylim(-0.5, len(msa) - 0.5)

            plt.tick_params(axis="both", which="both", labelsize=18)

            cbar = plt.colorbar()
            cbar.set_label("Sequence identity to query", fontsize=24, labelpad=24)
            cbar.ax.tick_params(labelsize=18)
            plt.xlabel("Positions", fontsize=24, labelpad=24)
            plt.ylabel("Sequences", fontsize=24, labelpad=36)
            plt.savefig(image_path)

        # Interactive plot of sequence coverage
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                z=final,
                colorscale="Rainbow_r",
                zmin=0,
                zmax=1,
                colorbar={"title": 'Your title'}
            )
        )
        # Add black line for sequence coverage depth
        fig.add_trace(
            go.Scatter(
                x=list(range(len(column_counts))),
                y=column_counts,
                mode="lines",
                line=dict(color="black", width=2),
                name="Coverage Depth",
            )
        )
        fig.update_layout(
            title=dict(text="Sequence coverage", x=0.5, xanchor="center"),
            xaxis_title="Positions", yaxis_title="Sequences",
        )

    if save_image:
        return fig, image_path
    else:
        return fig

def generate_pae_plot(pae_path, out_dir, name, save_image=True):
    """
    Generate a Plotly heatmap for Predicted Aligned Error (PAE) data.

    Args:
        pae (2D array): The PAE matrix.
    Returns:
        fig: A Plotly figure object of the PAE heatmap in green color scale
    """
    pae = np.genfromtxt(pae_path, delimiter="\t")
    max_pae = np.max(pae)
    fig = go.Figure()

    # Add heatmap
    fig.add_trace(
        go.Heatmap(
            z=pae,
            colorscale="Greens_r",
            zmin=0,
            zmax=max_pae,
        )
    )

    if save_image:
            image_path = f"{out_dir}/{name+('_' if name else '')}pae.png"
            fig.write_image(image_path, width=800, height=800)

    return fig
