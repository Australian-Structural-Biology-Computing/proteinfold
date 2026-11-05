#!/usr/bin/env python3
"""Capture PAE and pTM from the inference performed by the stock esm-fold CLI.

RUN_ESMFOLD copies this file to sitecustomize.py in its task directory and puts
that directory on PYTHONPATH only for the esm-fold command. Python imports
sitecustomize at startup, so the CLI and its single inference run as usual.
"""

import csv
import os
from pathlib import Path


pae_file = os.environ.get('ESMFOLD_PAE_FILE')
ptm_file = os.environ.get('ESMFOLD_PTM_FILE')

if pae_file or ptm_file:
    if not pae_file or not ptm_file:
        raise ValueError('Both ESMFOLD_PAE_FILE and ESMFOLD_PTM_FILE are required')

    from esm.esmfold.v1.esmfold import ESMFold

    _original_infer = ESMFold.infer

    def _infer_and_capture_metrics(self, sequences, *args, **kwargs):
        sequence_list = [sequences] if isinstance(sequences, str) else list(sequences)
        if len(sequence_list) != 1:
            raise ValueError('RUN_ESMFOLD expects one sequence per esm-fold invocation')

        output = _original_infer(self, sequences, *args, **kwargs)
        present = output['atom37_atom_exists'][0].any(dim=-1)
        pae = output['predicted_aligned_error'][0][present][:, present]
        expected_length = len(sequence_list[0].replace(':', ''))
        if tuple(pae.shape) != (expected_length, expected_length):
            raise ValueError(
                f'PAE shape {tuple(pae.shape)} does not match sequence length {expected_length}'
            )
        if not pae.isfinite().all():
            raise ValueError('ESMFold produced non-finite PAE values')

        ptm = output['ptm'][0]
        if not ptm.isfinite():
            raise ValueError('ESMFold produced a non-finite pTM value')

        with Path(pae_file).open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            for row in pae.detach().cpu():
                writer.writerow([f'{value:.4f}' for value in row.tolist()])
        with Path(ptm_file).open('w', newline='') as handle:
            csv.writer(handle, delimiter='\t').writerow(
                ['rank_0', f'{ptm.item():.6f}']
            )
        return output

    ESMFold.infer = _infer_and_capture_metrics
