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

User-facing outputs are largely consistent across modes.

<details markdown="1">
<summary>Common output patterns</summary>

- `<MODE>/top_ranked_structures/<SEQUENCE NAME>.{pdb,cif}` (format depends on the selected mode)
- `<MODE>/<SEQUENCE NAME>_colabfold2_metrics.jsonl.gz` (scalar metrics, one record per value)
- `<MODE>/<SEQUENCE NAME>_colabfold2_metrics.npz` (per-residue pLDDT and per-rank PAE matrices)
- `msa/<MSA TOOL>/<SEQUENCE NAME>.a3m` (raw MSA; see below)

</details>

### ColabFold2

ColabFold2 outputs are published directly under the selected public mode token, for example `alphafold3/<SEQUENCE NAME>/`, `boltz2/<SEQUENCE NAME>/`, or `opendde/<SEQUENCE NAME>/`. All three ESMFold2 model selections publish under `esmfold2/<SEQUENCE NAME>/`.

<details markdown="1">
<summary>Output files</summary>

- `<MODE>/` (the selected public mode token)
  - `top_ranked_structures/<SEQUENCE NAME>.cif` (`.pdb` for the `alphafold2` modes) (top-ranked structure per input)
  - `<first sample of the batch>_colabfold2_metrics.jsonl.gz` (scalar metrics for every sample in the batch)
  - `<first sample of the batch>_colabfold2_metrics.npz` (pLDDT vectors and PAE matrices for every sample in the batch)
  - `raw/` (published when `--save_intermediates` is enabled)

</details>

> [!NOTE]
> Select the ESMFold2 language-model variants with `--mode esmfold2 --esmfold2_model lm300m` or `--mode esmfold2 --esmfold2_model lm600m`. Confidence output depends on the checkpoint and runtime: pLDDT may contain finite scores or `n/a`, PAE may contain finite scores or `NaN`, and unavailable summary metrics remain empty.

### Packed metric stores (`*_colabfold2_metrics.jsonl.gz` / `.npz`)

All per-sequence metrics are written into two packed files per inference
batch instead of one small file per metric per sequence (which would exhaust
file-system inodes at genome scale). The batch is named after its first
sample; the stores hold every sample in the batch.

**`<first sample>_colabfold2_metrics.jsonl.gz`** - gzipped JSON Lines, one
record per scalar value:

```json
{"sample_id": "T1024", "rank": 1, "metric": "iptm", "chain_pair": null, "value": 0.712}
{"sample_id": "T1024", "rank": 1, "metric": "chainwise_ipsae", "chain_pair": "A:B", "value": 0.995}
{"sample_id": "T1024", "rank": 0, "metric": "msa_depth", "chain_pair": null, "value": 3}
```

`metric` is one of `ptm`, `iptm`, `ipsae`, `chainwise_iptm`,
`chainwise_ipsae` (the chainwise metrics carry a `chain_pair` label such as
`A:B`), `msa_depth` and `msa_length`. The `plddt_structure_only` record
(flags predictions without a confidence head) is described in the npz file
below.

**`<first sample>_colabfold2_metrics.npz`** - a compressed numpy archive
with one member per sequence and rank, keyed `<sample_id>|<rank>|<kind>`:

- `plddt`: 1-D float64 array of per-residue pLDDT values (rounded to 2
  decimals, in residue order).
- `pae`: 2-D float64 (n x n) PAE matrix per ranked model.

Read them with `numpy.load`; members decompress lazily. The stores are
consumed directly by the pipeline's own MultiQC summary, HTML reports and
ModelCIF assembly - no per-metric intermediate files are produced.

For quick inspection outside Python, the JSONL store is plain `zcat`/`grep`
text.

### Per-mode reports

<details markdown="1">
<summary>Output files</summary>

- `reports/`
  - `<SEQUENCE NAME>_<MODE>_report.html` (report per sequence)

</details>

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

- `msa/<MSA TOOL>/` (raw MSA search results such as `<SEQUENCE NAME>.a3m` and AlphaFold3-format JSON files; published by the `colabfold_search` tool)

### `--save_intermediates`

If `--save_intermediates` is enabled, extra raw intermediate files are published in mode-specific `raw/` directories.

Examples include:

- `<MODE>/raw/` for the selected ColabFold2 mode (scores JSONs, ranked structures and the a3m MSAs of the batch; flat, sample-prefixed names)

These raw outputs are intended for advanced debugging, reproducibility and method-specific downstream analyses. For detailed, canonical tool-specific native output specifications, see:

- [AlphaFold2](https://github.com/google-deepmind/alphafold?tab=readme-ov-file#alphafold-output)
- [AlphaFold3](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)
- [ColabFold](https://www.ebi.ac.uk/training/online/courses/alphafold/advanced-modeling-and-applications-of-predicted-protein-structures/customising-alphafold-structure-predictions/outputs-from-colabfold/)
