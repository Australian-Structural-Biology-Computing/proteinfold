# nf-core/proteinfold: Output

## Introduction

This document describes the user-facing output produced by the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and predicts protein structures using the following methods:

- [AlphaFold2](https://github.com/google-deepmind/alphafold)
- [AlphaFold3](https://github.com/google-deepmind/alphafold3)
- [ColabFold](https://github.com/sokrypton/ColabFold)

See main [README.md](https://github.com/nf-core/proteinfold/blob/master/README.md) for a condensed overview of the steps in the pipeline, and the bioinformatics tools used at each step.

The directories listed below will be created in the output directory after the pipeline has finished. All paths are relative to the top-level results directory.

Exact subdirectories depend on the selected mode. A run will typically contain top-level directories such as `top_ranked_structures/`, `metrics/`, `multiqc/`, `reports/`, and `pipeline_info/`.

### Prediction outputs (all modes)

The output directory stays flat: every result file carries a short mode tag
(`<MODE>`) and a short MSA-source tag (`<MSA>`) in its name instead of living
in per-mode directories, so multiple prediction modes and MSA sources can
write into the same tree without colliding. The mode shorthands are `af2`
(alphafold2), `af3` (alphafold3), `bz2` (boltz2), `chai1`, `esm2` (esmfold2),
`if2` (intellifold2), `odde` (opendde), `of3` (openfold3), `ob0` (openbind0),
`px2` (protenix2) and `rf3` (rosettafold3); all three ESMFold2 model
selections publish under the `esm2` tag. The MSA-source shorthands are
`mmseqs` (`colabfold_search`), `af2msa` (`alphafold2_msa`), `af3msa`
(`alphafold3_msa`), `server` (ColabFold MSA server) and `byo`
(sample-sheet-provided MSAs).

<details markdown="1">
<summary>Output files</summary>

- `top_ranked_structures/<SEQUENCE NAME>_<MODE>_<MSA>.{cif,pdb}` (`.pdb` for the `alphafold2` modes) (top-ranked structure per input, mode and MSA source)
- `metrics/<SEQUENCE NAME>_<MODE>_<MSA>_<METRIC>.tsv` with `<METRIC>` one of
  `plddt`, `confidence`, `ptm`, `iptm`, `ipsae`, `chainwise`, `msa`
  (per-prediction MSA depth summary)
- `metrics/<SEQUENCE NAME>_<MODE>_<MSA>_<RANK>_pae.tsv` (rank-1 only; published when the backend emits a PAE)
- `reports/<SEQUENCE NAME>_<MODE>_<MSA>_report.html` (report per sequence, mode and MSA source)
- `raw/` (published when `--save_intermediates` is enabled; file names embed the backend model type)

</details>

> [!NOTE]
> Select the ESMFold2 language-model variants with `--mode esmfold2 --esmfold2_model lm300m` or `--mode esmfold2 --esmfold2_model lm600m`. Confidence output depends on the checkpoint and runtime: pLDDT may contain finite scores or `n/a`, PAE may contain finite scores or `NaN`, and unavailable summary metrics remain empty.

### pLDDT (`<SEQUENCE NAME>_<MODE>_<MSA>_plddt.tsv`)

Confidence values per residue, rounded to 2 decimal places. Each ranked result gets its own column (for all-atom modules, atomic token confidences are processed to a naive mean value across the residue).

```
Positions	rank_0	rank_1	rank_2	rank_3	rank_4
0	83.58	85.27	88.41	86.22	84.91
1	97.99	97.81	97.39	97.49	97.32
2	98.22	98.42	98.16	97.88	97.81
3	98.06	98.15	97.94	97.56	97.4
4	98.67	98.56	98.3	98.38	98.29
5	98.81	98.77	98.62	98.61	98.54
6	98.79	98.74	98.57	98.59	98.52
...
```

### MSA (`<SEQUENCE NAME>_<MODE>_<MSA>_msa.tsv`)

The amino acid characters are converted to integers `0-19`, unknown as 20, **integer `21`** represents the gap character.

```
19	5	5	4	10	16	15	3	8	15	13	16	12	9	17	16	9	4	8	11	0	7	7	8	11	0	19	8	8	5	3
19	5	5	4	10	16	15	3	8	15	13	16	12	9	17	16	9	4	8	11	0	7	7	8	11	0	6	8	8	5	13
19	5	5	4	10	5	15	13	14	0	14	16	12	9	17	16	9	4	14	11	0	7	5	8	15	4	5	8	3	5	21
19	5	5	4	10	16	15	3	8	15	13	16	12	9	17	16	9	4	8	11	0	7	7	8	11	0	19	8	8	5	21
19	5	5	4	10	16	15	3	8	15	13	16	12	9	17	16	9	4	8	11	0	7	7	8	11	0	19	8	8	5	13
19	5	5	4	10	16	15	3	8	15	13	16	12	9	7	16	9	4	8	11	0	7	7	8	11	0	6	8	8	5	13
```

This allows easy sequence indentity calculation when processing as a `numpy` array.

### (i)pTM (`<SEQUENCE NAME>_<MODE>_<MSA>_[i]ptm.tsv`)

(i)pTM scores, rounded to 3 decimal places, listed by the rank number (currently unsorted - to reflect models and seeds where appropriate).

```
17  0.552
22  0.529
21  0.532
20  0.541
23  0.523
3 0.606
2 0.610
4 0.606
1 0.616
0 0.617
12  0.580
9 0.588
13  0.580
11  0.583
14  0.570
15  0.565
24  0.517
16  0.560
18  0.550
19  0.550
10  0.588
5 0.600
6 0.597
7 0.596
8 0.595
```

### chain-wise pTM/iPTM (`<SEQUENCE NAME>_<MODE>_<MSA>_chainwise.tsv`)

Chain-wise pTM and interface iPTM values, rounded to 4 decimal places, with chain-pair lettering as the row (`X:Y`) and rank number as the column. The on-diagonal `X:X` rows are the per-chain pTM; the off-diagonal `X:Y` rows are the pairwise interface iPTM. The values are program-native for every mode: the AlphaFold3-family backends write the chain-pair ipTM matrix and per-chain pTM of their own confidence head into the scores JSON (surfaced by the bundled `colabfold_batch_af3_chainwise.py` wrapper), and the AlphaFold2 backends write the same keys from their `--calc-extra-ptm` pass; both report one entry per unordered chain pair (the AF3 confidence head writes the same value into both orientations of its symmetric matrix). Ranks a mode does not score are `n/a`.

```
0	1	2
A:A	0.3810	0.3790	0.3850
A:B	0.2880	0.2750	0.2900
B:B	0.4000	0.3980	0.4050
```

In the HTML reports, chainwise iPTM and ipSAE are displayed as chain-by-chain matrices; the chainwise iPTM matrix carries the per-chain pTM on its diagonal (A:A). Modes that do not emit these metrics omit the corresponding report sections. Each report shows the top-ranked structure of each requested mode.

**ipSAE Attribution**

- The global ipSAE (and the per-chain-pair pDockQ/pDockQ2 scores) come from colabfold's own `get_interface_scores` implementation, which the AlphaFold2 backends run natively and the AlphaFold3-family backends receive through the bundled `colabfold_batch_af3_chainwise.py` wrapper. That implementation is a rework of the IPSAE project by the Dunbrack Lab (https://github.com/DunbrackLab/IPSAE/); consult the IPSAE repository for full details and citation information.

### chain-wise ipSAE (`<SEQUENCE NAME>_<MODE>_<MSA>_chainwise_ipsae.tsv`)

Per-chain-pair ipSAE values, rounded to 4 decimal places, with chain-pair lettering as the row (`X:Y`) and rank number as the column. These come straight from colabfold's native `get_interface_scores` computation, which the AlphaFold2 backends run for every complex and the AlphaFold3-family backends receive through the bundled `colabfold_batch_af3_chainwise.py` wrapper. ipSAE is asymmetric - the alignment runs from each chain against the other - so both orientations are kept:

```
0	1
A:B	0.7100	0.3050
B:A	0.6900	0.2980
```

In the HTML reports, the ipSAE matrix is displayed chain-by-chain for each ranked model alongside the chainwise iPTM matrix. Modes that do not emit this metric omit the report section.

### PAE (`<SEQUENCE NAME>_<MODE>_<MSA>_<RANK>_pae.tsv`)

Predicted alignment error of residues `j` aligned by residue `i`, rounded to 4 decimal places.
The row number gives you the index of residue `i` and the column value within the row gives the index of residue `j` for the 2D PAE matrix.

Each model prediction generates a separate file containing the rank number. Rank numbering follows the native convention of the underlying tool, so top-ranked models may appear as either `_0_pae.tsv` or `_1_pae.tsv` depending on the mode. Additional ranked results are stored within the `paes/` folder.

```
0.2500	1.5710	3.9037	6.2177	8.4471	11.4583	12.9679	15.1237	18.0263	18.3868	18.9381	20.5747	19.3314	20.1825	21.6145	23.2190
2.2177	0.2500	1.5559	4.0327	6.3151	7.6372	10.1969	11.3626	14.9366	16.1303	17.9119	19.1877	21.2715	20.9531	20.1760	19.4087
3.4270	1.5284	0.2500	2.1333	3.5351	5.1049	6.6521	8.2317	12.1379	13.7185	14.9523	16.6154	19.6988	21.7614	18.6592	17.9619
6.1051	5.4206	2.5987	0.2500	2.0724	5.1454	6.7492	9.5538	9.6285	12.3868	13.8527	16.3586	17.2605	20.6381	19.9987	19.3295
7.3512	6.4947	5.5435	2.6740	0.2500	1.7561	4.9041	6.3923	8.9735	8.9272	12.3419	14.6005	15.9820	17.6358	20.5190	19.1028
7.4734	7.0899	5.8128	5.7512	2.0439	0.2500	1.8352	5.1064	6.4225	9.2098	10.5136	12.9404	14.3152	16.8122	18.6336	17.7382
```

#### Example report plots

The report shows the top-ranked structure in a Molstar viewer beside the Predicted Aligned Error (PAE) matrix, with the MSA sequence-coverage plot underneath. The views are cross-linked: drag-selecting a block in the PAE matrix or a position range in the MSA coverage plot highlights the corresponding residues in the structure, and clicking a residue in the structure highlights it in the plots. The PAE heatmap is resampled with a proper area average (each plot pixel averages every token whose lane overlaps it), so blocks and diagonal stripes stay faithful when the matrix is drawn smaller than the token count. Hovering the PAE reports the expected position error in Ångströms for the residue pair; hovering the MSA reports sequence identity to the query and column coverage. A sequence panel shows each chain's residues as one-letter codes in a numbered grid whose row length adapts to the panel width in steps of ten (rows break at chain boundaries; ligands keep their residue names and are coloured by their chain id like every other chain), with the chain identity shown only by a coloured chip on a heading line above each chain (a Hide/Show toggle in the section header collapses the panel); clicking a residue or dragging a range there highlights it in the structure, the PAE matrix and the MSA coverage plot (Shift adds to the selection instead of replacing it, and removes a click or range that is already fully selected), and selections made in any other view highlight the corresponding residues in the sequence with a continuous band under the covered letters. Black separator lines mark the chain boundaries on both plots, and the MSA rows are ordered with the query (and the most similar sequences) at the top, drawn in the same rainbow colour scale as the previous static coverage plot. The MSA coverage plot embeds every sequence (no row subsampling) and resamples the full matrix into the plot area. In the ipSAE/ipTM colourings a click selects whole chains instead of single residues: a diagonal block selects that chain and an off-diagonal block selects both chains of the pair. Holding Shift while clicking or dragging in the PAE matrix or the MSA coverage plot adds the new selection to the current one instead of replacing it. Those colourings are click-only (dragging selects nothing), the selection itself persists across colouring switches (the structure, MSA and sequence keep their highlights), but residue-level bands are drawn on the matrix only in the PAE colouring - switching back to it re-draws the current selection with the same row/column bands it had before. Residues are coloured by pLDDT (from the B-factor column); a `pLDDT | Chain` toggle above the viewer (selected tabs are marked with a dark pill) switches to chain colouring where every atom of a chain stays within that chain's colour family — carbons take the chain colour itself (the same shared palette as the chain legend chips and the matrix axis bands, so each chain reads as the same colour across views) while every other atom keeps its default CPK element colour; ligand chains are treated identically, and the click-to-focus visual of a selected residue follows the active colouring as well.

Additional PAE features:

- **Ligand lanes:** non-protein tokens (ligands, cofactors, ions) are grouped into fixed-width lanes on both PAE axes, coloured by their chain id like the polymer chains; each chain of a multi-chain ligand prediction gets its own lane. The matrix works in the prediction's native token space: per-token chain ids from the prediction output (written into the PAE TSV by the metrics extraction) place chain boundaries and ligand lanes correctly even when ligands are scored one token per atom, and every token is mapped back onto its structure residue so cross-highlighting with the structure, MSA and sequence keeps working. Only a matrix without any token information falls back to the unlabelled mismatch warning. Hovering a ligand lane names the ligand and its PAE to the selected residue. The matrix is drawn in the standard PAE orientation: the first token sits in the top-left corner, so chain bands read in the same order on both axes.
- **Off-diagonal selection:** dragging a rectangle that spans two different chain blocks selects both the aligned (row) block and the scored (column) block — drawn as a horizontal band across the selected rows crossing a vertical band down the selected columns; a single-cell click selects the residue pair. Selections made from the structure or the MSA draw both bands through each selected residue.- **Chainwise score overlay:** when chainwise ipSAE or ipTM values are available, a `PAE | ipTM | ipSAE` toggle above the matrix switches the block colouring from the PAE heatmap to the interface score of each chain pair (same-chain blocks are grey, unreported pairs near-white), with the score value printed on each reported block - in a label colour chosen by the block fill's luminance so it stays readable across the whole score ramp - and only when the block is large enough to hold it: small blocks such as the fixed-width ligand lanes stay clean and report their score on hover. The tooltips switch accordingly and the fallback chainwise tables only appear when no PAE matrix is available.

##### Sequence coverage

![Sequence coverage](images/sequence_coverage_proteinfold-v2.png?raw=true "Example sequence coverage plot")

##### Predicted Aligned Error (PAE)

![PAE](images/pae_proteinfold-v2.png?raw=true "Example PAE plot")

> The example figures above predate the current report layout, which shows the
> Molstar structure viewer beside the PAE matrix with the pLDDT colour legend
> instead of a separate per-residue pLDDT plot.

### Per-mode reports

<details markdown="1">
<summary>Output files</summary>

- `reports/`
  - `<SEQUENCE NAME>_<MODE>_<MSA>_report.html` (report per sequence, mode and MSA source; both are shorthand tags)

</details>

### Foldseek structural similarity search

### MultiQC report

<details markdown="1">
<summary>Output files</summary>

- `multiqc`
  - `*_multiqc_report.html`: Standalone HTML report(s) that can be viewed in your web browser.
  - `*_multiqc_report_data/`: Parsed report data for each corresponding MultiQC report.

</details>

[MultiQC](https://multiqc.info/docs/) is a visualisation tool that generates HTML report(s) summarising samples in your project. Most QC results are visualised in the report and further statistics are available within each corresponding `*_multiqc_report_data/` directory.

Results generated by MultiQC collate QC metrics from the selected structure-prediction mode(s), and the software versions for traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### Additional intermediate outputs

Depending on the selected options, additional top-level directories may be present, for example:

- `msa/` (raw MSA results; file names carry an MSA-tool shorthand — `mmseqs` for `colabfold_search`, `af2msa`/`af3msa` for the data-pipeline tools, e.g. `<SEQUENCE NAME>_mmseqs.a3m` and `<SEQUENCE NAME>_mmseqs.json`)

### `--save_intermediates`

If `--save_intermediates` is enabled, extra raw intermediate files are published under a flat `raw/` directory.

Examples include:

- `raw/` for every selected ColabFold2 mode (the file names embed the backend model type, e.g. `<SEQUENCE NAME>_rank_001_alphafold3_seed_000_sample_0.cif`)

These raw outputs are intended for advanced debugging, reproducibility and method-specific downstream analyses. For detailed, canonical tool-specific native output specifications, see:

- [AlphaFold2](https://github.com/google-deepmind/alphafold?tab=readme-ov-file#alphafold-output)
- [AlphaFold3](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)
- [ColabFold](https://www.ebi.ac.uk/training/online/courses/alphafold/advanced-modeling-and-applications-of-predicted-protein-structures/customising-alphafold-structure-predictions/outputs-from-colabfold/)
