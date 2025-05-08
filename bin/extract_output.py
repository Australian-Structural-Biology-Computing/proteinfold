#!/usr/bin/env python
import pickle
import os
import argparse
import json
import torch
import numpy as np
import csv
from utils import plddt_from_struct_b_factor

# Mapping of characters to integers
AA_to_int = {
    "A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9,
    "M": 10, "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19,
    "-": 21, ".": 21
}

def extract_structs_plddt_to_tsv(id, struct_files):
    """
    Write out a tsv file for reading by MultiQC in nf-core/proteinfold
    Uses utils function with BioPython PDB package to extract residue pLDDT values from the b-factor column.
    """

    # Set up headers
    with open(id + '_plddt_mqc.tsv', 'w') as multiqc_tsv:
        writer = csv.writer(multiqc_tsv, delimiter='\t')
        rank_names = []
        for i in range(len(struct_files)):
            rank_names.append(f"rank_{i}")
        writer.writerow(["Positions"] + rank_names)

    res_counts = []
    plddt_cols = []

    for struct_file in struct_files:
        plddts = plddt_from_struct_b_factor(struct_file)
        res_counts.append(len(plddts))
        plddt_cols.append(plddts)

    if not all(x == res_counts[0] for x in res_counts):
        print("Not all structures have the same number of residues!")
        return
    else:
        res_id_col = list(range(len(plddts)))

    # Check the pLDDT cols are the same size before combining
    if not (len(set(len(col) for col in plddt_cols)) == 1):
        print("Not all pLDDT columns have the same number of values!")

    plddt_rows = zip(res_id_col, *plddt_cols) #combine lists column-wise to make rows

    with open(id + '_plddt_mqc.tsv', 'a') as multiqc_tsv:
        writer = csv.writer(multiqc_tsv, delimiter='\t')
        for plddt_row in plddt_rows:
            writer.writerow(plddt_row)


def read_pkl(id, pkl_files):
    """
    Adapted from the Galaxy AlphaFold tool (https://github.com/usegalaxy-au/tools-au/blob/de94df520c8dc7b8652aedb92e90f6ebb312f95f/tools/alphafold/scripts/outputs.py), originally authored by @neoformit and @graceahall and funded by Australian Biocommons and QCIF Australia.
    """

    for pkl_file in pkl_files:
        print(pkl_file)
        data = pickle.load(open(pkl_file, "rb"))
        if pkl_file.endswith("final_features.pkl"): #HelixFold3
            with open(f"{id}_msa.tsv", "w") as out_f:
                for val in data["feat"]["msa"]:
                    out_f.write("\t".join([str(x) for x in val]) + "\n")  # TODO: take this out as a line
        elif pkl_file.endswith("features.pkl"):  #AlphaFold2.3
                with open(f"{id}_msa.tsv", "w") as out_f:
                    for val in data["msa"]:
                        out_f.write("\t".join([str(x) for x in val]) + "\n")
        else:  #AlphaFold2.3 non-summary. TODO: Need to either read in ranking_debug.json to get the ranking order, or do it later in the workflow.
            model_id = (
                os.path.basename(pkl_file)
                .replace("result_model_", "")
                .replace(".pkl", "")
            )
            with open(f"{id}_lddt_{model_id}.tsv", "w") as out_f:
                out_f.write("\t".join([str(x) for x in data["plddt"]]) + "\n")

            # Full credit to Cam Hyde @ QCIF/Galaxy who I lifted this from: https://github.com/usegalaxy-au/tools-au/blob/de94df520c8dc7b8652aedb92e90f6ebb312f95f/tools/alphafold/scripts/outputs.py
            if 'predicted_aligned_error' not in data.keys():
                print(f"No PAE output in {pkl_file}, it was likely a monomer calculation")
            else:
                with open(f"{id}_pae_{model_id}.tsv", "w") as out_f:
                    PAE = data["predicted_aligned_error"]
                    for row in PAE:
                        rounded_row = [f"{num:.4f}" for num in row]
                        out_f.write('\t'.join(rounded_row) + '\n')

def a3m_to_int(a3m_file):  # For the RosettaFold-All-Atom .a3m. Written with GitHub Copilot
    """
    Convert an A3M MSA representation into an integer representation (0-21).
    """
    # Tom Litfin gave me format run-down, the lowercase are just insertions

    with open(a3m_file, "r") as f:
        msa = f.read()

    # Convert each sequence in the MSA
    int_sequences = []
    for idx, line in enumerate(msa.splitlines()):
        if idx == 0 and not line.startswith(">"):  # If there's an additional header (non-FASTA) skip it. E.g ColabFold
            continue

        if not line.startswith(">"):  # Ignore header lines
            filtered_line = ''.join(char for char in line if not char.islower()) # Remove inserts (lower-case chars) in a3m
            int_sequence = [AA_to_int.get(char.upper(), 20) for char in filtered_line]
            int_sequences.append(int_sequence)

    int_sequences_array = np.array(int_sequences, dtype=object)


    return int_sequences_array

def read_a3m(id, a3m_files):
    for a3m_file in a3m_files:
        int_seqs = a3m_to_int(a3m_file)
        with open(f"{id}_msa.tsv", "w") as out_f:
                for row in int_seqs:
                    out_f.write("\t".join(map(str, row)) + "\n")

def read_npz(id, npz_files):
   for npz_file in npz_files:
       if npz_file.split('/')[-1].startswith('pae') and npz_file.endswith('.npz'): #Boltz PAE files if --write_full_pae is used
            data = np.load(npz_file)  # using a with open()
            PAE = data['pae']
            with open(f"{id}_pae.tsv", "w") as out_f:
                for row in PAE:
                    rounded_row = [f"{num:.4f}" for num in row]
                    out_f.write('\t'.join(rounded_row) + '\n')

#MSA data is here: https://github.com/jwohlwend/boltz/blob/1f7acb18f279858bc292a8a0f9fbb5d96d6491f1/src/boltz/data/types.py#L298-L315  it looks like with  ("res_type", np.dtype("i1")) it's an undending list,  and you need "sequences" to get start and end indices

#       msa = data["residues"]  # TODO - dump these intio to a row seperate MSA like the others
#       msa_2 = data["sequences"]  #sequences just seems to list length
#       print(msa)
#       print(msa_2)

#       with open(f"{id}_msa.tsv", "w") as out_f:
#           for val in data["sequences"]:
#               out_f.write("\t".join([str(x) for x in val]) + "\n") # TODO - this has a single line for each entry, fix it up


def read_json(id, json_files):
    for json_file in json_files:
        if json_file.endswith("_data.json"): #AF3 output with MSA info
            with open(json_file, 'r') as f:
                data = json.load(f)
                unpaired_MSAs = data['sequences'][0]['protein']['unpairedMsa']

                msa_lines = [line for line in unpaired_MSAs.split("\n") if not line.startswith(">") and line.strip()]

                int_seqs = [[AA_to_int.get(residue, 20) for residue in line] for line in msa_lines]

                with open(f"{id}_msa.tsv", "w") as out_f:
                        for row in int_seqs:
                            out_f.write("\t".join(map(str, row)) + "\n")

        elif json_file.endswith("_confidences.json") or json_file.endswith('all_results.json'): #AF3 output with PAE info, or HF3 PAE data. TODO: Need to make sure the workflow points to [protein]/[protein]_rank1/all_results.json
            with open(json_file, 'r') as f:
                data = json.load(f)
                PAE = data['pae']

                with open(f"{id}_pae.tsv", "w") as out_f:
                    for row in PAE:
                        out_f.write('\t'.join([str(round(x,4)) for x in row]) + '\n')  #tsv since the other metrics are .tsv in proteinfold

        elif json_file.endswith("predicted_aligned_error_v1.json"): #ColabFold file
            with open(json_file, 'r') as f:
                data = json.load(f)
                PAE = data['predicted_aligned_error']

                with open(f"{id}_pae.tsv", "w") as out_f:
                    for row in PAE:
                        rounded_row = [f"{num:.4f}" for num in row]
                        out_f.write('\t'.join(rounded_row) + '\n')

def read_pt(id, pt_files):
    for pt_file in pt_files:
            with open(pt_file, 'rb') as f:   # TODO: point to [protein]_aux.pt
                data = torch.load(f, map_location="cpu")
                PAE = data['pae']

                with open(f"{id}_pae.tsv", "w") as out_f:
                    for tensor in PAE.tolist():
                        for row in tensor:
                            rounded_row = [f"{num:.4f}" for num in row]
                            out_f.write('\t'.join(rounded_row) + '\n')

parser = argparse.ArgumentParser()
parser.add_argument("--pkls", dest="pkls", required=False, nargs="+") # For reading both HelixFold3 and AlphaFold2 MSA formats
parser.add_argument("--npzs", dest="npzs", required=False, nargs="+") # For reading the Boltz-1 PAE formats. TODO: Boltz-1 MSA not implemented (go straight to .a3m file), implement
parser.add_argument("--a3ms", dest="a3ms", required=False, nargs="+") # For reading the RosettaFold-All-Atom, ColabFold, and Boltz-1 MSA formats
parser.add_argument("--jsons", dest="jsons", required=False, nargs="+") # For reading the AF3 MSA & PAE, HF3 PAE
parser.add_argument("--pts", dest="pts", required=False, nargs="+") # For read RFAA pytorch model to get PAE data
parser.add_argument("--structs", dest="structs", required=False, nargs="+")
parser.add_argument("--name", default="untitled", dest="name") # might need a --name $meta.id
parser.add_argument("--output_dir", default=".", dest="output_dir")
args = parser.parse_args()

if args.pkls is not None:
    read_pkl(args.name, args.pkls)
if args.a3ms is not None:
    read_a3m(args.name, args.a3ms)
if args.npzs is not None:
    read_npz(args.name, args.npzs)
if args.jsons is not None:
    read_json(args.name, args.jsons)
if args.pts is not None:
    read_pt(args.name, args.pts)
if args.structs is not None:
    extract_structs_plddt_to_tsv(args.name, args.structs)
