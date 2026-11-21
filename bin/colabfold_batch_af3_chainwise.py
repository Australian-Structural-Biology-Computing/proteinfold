#!/usr/bin/env python3
"""colabfold_batch wrapper that surfaces the AF3-family chainwise metrics.

The af3-preview colabfold backend computes chain-level confidence data
during inference -- ``chain_pair_iptm`` (pairwise interface ipTM) and
``iptm_ichain`` (per-chain pTM) are present in the inference result
metadata -- but its ``scores_of`` leaves them out of the scores JSON
(https://github.com/sokrypton/ColabFold/tree/af3-preview). Until that is
surfaced upstream, this wrapper patches
``colabfold.alphafold3.predict.scores_of`` to add

* ``pairwise_iptm``: one entry per unordered chain pair, the
  convention of the AlphaFold2 ``--calc-extra-ptm`` path (the AF3
  confidence head writes the same value into both orientations of its
  symmetric matrix),
* ``per_chain_ptm``: the per-chain pTM from ``iptm_ichain``,
* ``ipsae``/``pdockq``/``pdockq2``: colabfold's native per-chain-pair
  interface scores (``colabfold.alphafold.ipsae.get_interface_scores``,
  the same machinery the AlphaFold2 backend already runs), fed with the
  token-level PAE, the per-token chain segmentation and the CA/CB
  geometry recovered from the predicted structure,

under the keys the AlphaFold2 path uses, then hands off to colabfold's
own ``main`` so the CLI behaves identically.

Models whose inference results carry no chain summary in their metadata
(e.g. confidence-less models) are left un-augmented.
"""
import numpy as np

import colabfold.alphafold3.predict as af3_predict

_original_scores_of = af3_predict.scores_of


def _chain_labels(inference_result, summary_chain_ids):
    """Map the summary's chain order onto the predicted structure's chain ids.

    The structure chain ids are what the published structures and the
    interface scores use, so preferring them keeps every chainwise
    representation consistent. Falls back to the summary's own ids, then
    to A/B/C letters, when the structure's ids do not line up one-to-one
    with the summary order (e.g. duplicated chain ids).
    """
    structure_chain_ids = []
    for chain_id in np.asarray(inference_result.predicted_structure.chain_id):
        chain_id = str(chain_id)
        if chain_id not in structure_chain_ids:
            structure_chain_ids.append(chain_id)
    if len(structure_chain_ids) == len(summary_chain_ids):
        return structure_chain_ids
    if len(set(summary_chain_ids)) == len(summary_chain_ids):
        return [str(chain_id) for chain_id in summary_chain_ids]
    import string
    return [string.ascii_uppercase[i % len(string.ascii_uppercase)] for i in range(len(summary_chain_ids))]


def _af3_interface_scores(inference_result):
    """Run colabfold's native interface scores over the AF3 tokens.

    ``get_interface_scores`` consumes AlphaFold2-style arrays: a
    token-level PAE, a per-token chain id, and per-token CA/CB positions
    with a placement mask. The AF3 tokens carry their chain and residue
    ids in the confidence metadata, so the geometry is recovered from the
    predicted structure: a token contributes to the contact statistics
    when its residue owns a CA (or CB) atom; tokens without one (ligands,
    nucleic acids) stay in the chain segmentation for the PAE-based
    ipSAE but are unplaced for the distance-based pDockQ scores.
    """
    from alphafold.common import residue_constants
    from alphafold3.model import confidence_types
    from colabfold.alphafold import ipsae

    full = confidence_types.StructureConfidenceFull.from_inference_result(inference_result)
    struc = inference_result.predicted_structure

    ca_atom = residue_constants.atom_order["CA"]
    cb_atom = residue_constants.atom_order["CB"]
    num_atom_types = residue_constants.atom_type_num

    ca_map, cb_map = {}, {}
    atom_names = [str(name) for name in np.asarray(struc.atom_name)]
    atom_chains = [str(chain) for chain in np.asarray(struc.chain_id)]
    atom_res_ids = [int(res_id) for res_id in np.asarray(struc.res_id)]
    for index, name in enumerate(atom_names):
        key = (atom_chains[index], atom_res_ids[index])
        if name == "CA":
            ca_map[key] = index
        elif name == "CB":
            cb_map[key] = index

    token_chain_ids = list(full.token_chain_ids)
    token_res_ids = [int(res_id) for res_id in full.token_res_ids]
    num_tokens = len(token_chain_ids)
    chain_order = list(dict.fromkeys(token_chain_ids))
    chain_to_int = {chain: index + 1 for index, chain in enumerate(chain_order)}
    asym_id = np.array([chain_to_int[chain] for chain in token_chain_ids])

    atom_positions = np.zeros((num_tokens, num_atom_types), dtype=np.float64)
    atom_mask = np.zeros((num_tokens, num_atom_types), dtype=np.float64)
    plddt = np.zeros(num_tokens, dtype=np.float64)
    atom_x = np.asarray(struc.atom_x, dtype=np.float64)
    atom_y = np.asarray(struc.atom_y, dtype=np.float64)
    atom_z = np.asarray(struc.atom_z, dtype=np.float64)
    atom_b_factor = np.asarray(struc.atom_b_factor, dtype=np.float64)
    for token in range(num_tokens):
        key = (token_chain_ids[token], token_res_ids[token])
        ca = ca_map.get(key)
        if ca is None:
            continue
        atom_positions[token, ca_atom] = (atom_x[ca], atom_y[ca], atom_z[ca])
        atom_mask[token, ca_atom] = 1.0
        cb = cb_map.get(key)
        if cb is not None:
            atom_positions[token, cb_atom] = (atom_x[cb], atom_y[cb], atom_z[cb])
            atom_mask[token, cb_atom] = 1.0
        plddt[token] = atom_b_factor[ca]

    return ipsae.get_interface_scores(
        pae=full.pae,
        plddt=plddt,
        asym_id=asym_id,
        atom_positions=atom_positions,
        atom_mask=atom_mask,
    )


def scores_of_with_chainwise(inference_result, ranking_score):
    scores = _original_scores_of(inference_result, ranking_score)
    if scores.get("failed"):
        return scores
    try:
        from alphafold3.model import confidence_types

        summary = confidence_types.StructureConfidenceSummary.from_inference_result(inference_result)
        labels = _chain_labels(inference_result, summary.chain_ids)
        pair_iptm = np.asarray(summary.chain_pair_iptm, dtype=float)
        chain_ptm = np.asarray(summary.chain_ptm, dtype=float)
        num_chains = len(labels)
        if num_chains > 1 and pair_iptm.shape == (num_chains, num_chains) and chain_ptm.shape == (num_chains,):
            # One entry per unordered pair, as the AlphaFold2 path
            # reports them; the AF3 matrix is symmetric by construction,
            # so the max over the two halves is the shared value.
            scores["pairwise_iptm"] = {
                f"{labels[i]}-{labels[j]}": round(max(float(pair_iptm[i, j]), float(pair_iptm[j, i])), 3)
                for i in range(num_chains) for j in range(i + 1, num_chains)
            }
            scores["per_chain_ptm"] = {
                labels[i]: round(float(chain_ptm[i]), 3) for i in range(num_chains)
            }
        scores.update(_af3_interface_scores(inference_result))
    except (ImportError, KeyError, TypeError, ValueError):
        # No chain summary in this model's inference results: keep the
        # un-augmented scores rather than failing the prediction.
        pass
    return scores


af3_predict.scores_of = scores_of_with_chainwise


def main():
    from colabfold.batch import main as batch_main

    batch_main()


if __name__ == "__main__":
    main()
