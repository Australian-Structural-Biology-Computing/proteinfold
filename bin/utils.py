from Bio import PDB
import numpy as np

def plddt_from_struct_b_factor(struct_file):
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
