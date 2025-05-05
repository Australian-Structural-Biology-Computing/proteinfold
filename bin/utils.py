from collections import OrderedDict
import plotly.graph_objects as go
from Bio import PDB
import numpy as np

def pLDDT_from_struct_b_factor(struct_file):
    """
    Uses the BioPython PDB package to extract residue pLDDT values from the b-factor column. Iterates over PDB objects rather than processes raw file
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
    res_pLDDTs = []
    pLDDT_tot = 0

    for model in structure:
        for chain in model:
            chain_res_list = chain.get_unpacked_list()
            res_list.extend(chain_res_list)
            for residue in chain:
                atom_list = residue.get_unpacked_list()
                atom_pLDDT_tot = 0
                for atom in residue:  # ESMFold and others have separate atom-wise values, so doing atom-wise to cover that and residue-wise
                    atom_pLDDT = atom.get_bfactor()
                    atom_pLDDT_tot += atom_pLDDT

                res_pLDDT = float(atom_pLDDT_tot / len(atom_list))

                if (res_pLDDT < 1):  # RFAA the multiplication of mean isn't failing. Anyway covering to a [0,100] range for any structure file1
                    res_pLDDT *= 100

                res_pLDDTs.append(res_pLDDT)
                pLDDT_tot += res_pLDDT

    res_pLDDTs = np.array(res_pLDDTs)
    res_pLDDTs = np.round(res_pLDDTs, 2)

    return res_pLDDTs



def generate_plddt_plot(structures):
    """
    Generate a Plotly figure for predicted LDDT per position for given structures.

    Args:
        structures (list): List of structure file paths.

    Returns:
        go.Figure: Plotly figure object with LDDT data.
    """
    plddt_per_struct = OrderedDict()

    # Extract pLDDT values for each structure
    for struct in structures:
        plddt_per_struct[struct] = pLDDT_from_struct_b_factor(struct)

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
        title=dict(text="Predicted LDDT per position", x=0.5, xanchor="center"),
        xaxis=dict(
            title="Positions", showline=True, linecolor="black", gridcolor="WhiteSmoke"
        ),
        yaxis=dict(
            title="Predicted LDDT",
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
