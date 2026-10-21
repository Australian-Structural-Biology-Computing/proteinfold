#!/usr/bin/env python3
"""Run one ESMFold inference and save its structure, PAE, and pTM together."""

import argparse
import csv
from pathlib import Path

import esm
import torch
from esm.data import read_fasta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--fasta', type=Path, required=True)
    parser.add_argument('-o', '--pdb', type=Path, required=True)
    parser.add_argument('-m', '--model-dir', type=Path, required=True)
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--num-recycles', type=int, default=4)
    parser.add_argument('--chunk-size', type=int, default=None)
    parser.add_argument('--cpu-only', action='store_true')
    parser.add_argument('--cpu-offload', action='store_true')
    args = parser.parse_args()

    records = list(read_fasta(args.fasta))
    if len(records) != 1:
        raise ValueError(f'Expected one sequence in {args.fasta}, found {len(records)}')
    if args.cpu_only and args.cpu_offload:
        parser.error('--cpu-only and --cpu-offload cannot be used together')

    torch.hub.set_dir(str(args.model_dir))
    model = esm.pretrained.esmfold_v1().eval()
    model.set_chunk_size(args.chunk_size)
    if args.cpu_only:
        model.esm.float()
        model.cpu()
    elif args.cpu_offload:
        from esm.scripts.fold import init_model_on_gpu_with_cpu_offloading
        model = init_model_on_gpu_with_cpu_offloading(model)
    else:
        model.cuda()

    output = model.infer(records[0][1], num_recycles=args.num_recycles)
    output = {key: value.cpu() for key, value in output.items()}
    # ESMFold inserts an internal glycine linker for ':'-separated chains.
    # The PDB omits those residues, so remove the same rows/columns from PAE.
    present = output['atom37_atom_exists'][0].any(dim=-1)
    pae = output['predicted_aligned_error'][0][present][:, present]
    sequence_length = len(records[0][1].replace(':', ''))
    if tuple(pae.shape) != (sequence_length, sequence_length):
        raise ValueError(f'PAE shape {tuple(pae.shape)} does not match sequence length {sequence_length}')
    if not torch.isfinite(pae).all():
        raise ValueError('ESMFold produced non-finite PAE values')

    args.pdb.mkdir(parents=True, exist_ok=True)
    (args.pdb / f'{args.output_prefix}_esmfold.pdb').write_text(model.output_to_pdb(output)[0])
    with (args.pdb / f'{args.output_prefix}_pae.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for row in pae:
            writer.writerow([f'{value:.4f}' for value in row.tolist()])
    with (args.pdb / f'{args.output_prefix}_ptm.tsv').open('w', newline='') as handle:
        csv.writer(handle, delimiter='\t').writerow(['rank_0', f"{output['ptm'][0].item():.6f}"])


if __name__ == '__main__':
    main()
