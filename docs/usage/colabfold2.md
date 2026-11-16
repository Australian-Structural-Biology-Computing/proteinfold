---
title: ColabFold2
weight: 25
---

# ColabFold2

| Mode                                                 | Protein | RNA | Small-molecule | PTM | Constraints | pLM | MSA server | Split MSA |
| :--------------------------------------------------- | :-----: | :-: | :------------: | :-: | :---------: | :-: | :--------: | :-------: |
| [ColabFold2](https://github.com/sokrypton/ColabFold) |   ✅    | ✅  |       ✅       | ✅  |     ❌      | ✅  |     ✅     |    ✅     |

> [!WARNING]
> ColabFold2 support is experimental and its upstream interfaces may change. Which entities a backend actually supports is decided by that backend.

## General Usage

A ColabFold2 backend is selected with `--mode` (the token is the upstream backend name). Only one ColabFold2 backend can be run per invocation.

| Mode token     | Upstream model type                                                                                     |
| -------------- | ------------------------------------------------------------------------------------------------------- |
| `alphafold3`   | `alphafold3`                                                                                            |
| `boltz2`       | `boltz2`                                                                                                |
| `esmfold2`     | `esmfold2*`                                                                                             |
| `protenix2`    | `protenix2`                                                                                             |
| `chai1`        | `chai1`                                                                                                 |
| `intellifold2` | `intellifold2`                                                                                          |
| `opendde`      | `opendde`                                                                                               |
| `openfold3`    | `openfold3`                                                                                             |
| `openbind0`    | `openbind0`                                                                                             |
| `rosettafold3` | `rosettafold3`                                                                                          |
| `alphafold2`   | `auto` (colabfold_batch picks `alphafold2_ptm` for monomers and `alphafold2_multimer_v3` for complexes) |

```console
nextflow run nf-core/proteinfold \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --mode alphafold3 \
    --colabfold_db <null (default) | PATH> \
    --use_gpu \
    -profile <docker/singularity/podman/shifter/charliecloud/institute>
```

MSA generation is selected with `--msa_tool` (`colabfold_search`, `alphafold2_msa` or `alphafold3_msa`) or delegated to the ColabFold MSA server with `--use_msa_server`. The default `colabfold_search` runs the local MMseqs2 protocol against the databases in `--colabfold_db`; the `alphafold2_msa` and `alphafold3_msa` tools run the respective upstream data pipelines and convert their MSAs into the AlphaFold3-format input JSON used by the backends (see the [usage documentation](https://nf-co.re/proteinfold/usage)). `--db`/`--colabfold_db` is a database-only location; it is not expected to contain a `colabfold2_models` child when it is an S3 URL. If no explicit local `--colabfold2_params_path` is supplied, the pipeline reuses complete model weights under a local `<colabfold_db>/colabfold2_models` directory or runs the pinned upstream ColabFold2 downloader in a separate Nextflow task and stages the result in the task work directory. Outdated Chai1 or Boltz2 completion markers trigger a fresh download instead of reusing old blobs. For RosettaFold3 fp32, the pipeline also checks the checkpoint size pinned to this image: a stale implicit cache triggers a task-local download, while a stale explicit `--colabfold2_params_path` fails with a request to provide refreshed weights; neither case changes the source cache. This requires outbound network access and several GB of task storage; no model weights are downloaded in stub mode. AlphaFold3 JSON inputs are passed through to carry rich entities such as RNA, ligands, PTMs, templates and user CCD definitions: with offline MSA generation the pipeline extracts their protein queries, runs the search and merges the resulting MSAs back into the original JSON, while `--use_msa_server` (with an optional private `--msa_server_url`) passes the file directly to `colabfold_batch`.

The `alphafold2` modes run the classic AlphaFold2 models through `colabfold_batch`; they emit PDB instead of mmCIF structures and consume ColabFold a3m queries rather than the AlphaFold3 input JSON, so the pipeline converts the MSA JSON produced by the selected `--msa_tool` before inference. These models do not provide int8 weights; `--colabfold2_weights_precision int8` is rejected. With `--use_msa_server`, the server performs the search and the queries are passed as FASTA (AlphaFold3 JSON inputs are reduced to their protein chains first).

> [!NOTE]
> Local MSA generation occurs in a separate module to model inference, so inference can be re-run without repeating the MSA search.

## File Structure

If `--colabfold2_params_path` is not set and no weights cache is available under a local `--db`/`--colabfold_db`, the pipeline downloads the required weights in a separate Nextflow process. The weights root has this structure:

<details markdown="1">
<summary>Directory structure</summary>

```
<colabfold2_params_path>/params/af3/
├── libcifpp/
│   └── components.cif
├── esmc_300m|esmc_600m/
│   └── <variant>.bin.zst
└── <model>[-int8]/
    ├── <model>.bin.zst
    ├── <model>.lm.npz
    └── download_<model>_<precision>_finished.txt
```

The `alphafold2` modes use ColabFold's own AlphaFold2 parameter layout instead (`<colabfold2_params_path>/params/` with the `alphafold_params_*.npz` files written by `colabfold.download.download_alphafold_params`).

Chai1 and Boltz2 use revision-qualified completion markers instead: `download_chai1_<precision>_r3_finished.txt` and `download_boltz2_<precision>_r1_finished.txt`.

</details>

`<model>` is the selected upstream model, suffixed with `-int8` when `--colabfold2_weights_precision int8` is set; the `<model>.lm.npz` checkpoint and matching `esmc_300m` or `esmc_600m` trunk are only required by the ESMFold2 language-model selections.

## Additional Arguments

See the [ColabFold2](https://github.com/sokrypton/ColabFold) documentation for a full description of additional arguments. The arguments supported by the proteinfold workflow are described briefly below:

| Parameter                        | Default | Description                                                                                                                                                                                                                                  |
| -------------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--colabfold2_weights_precision` | `fp32`  | Checkpoint precision: `fp32` or `int8`. Official AlphaFold3 weights support `fp32` only, and the AlphaFold2 models provide no int8 weights.                                                                                                  |
| `--esmfold2_model`               | `msa`   | ESMFold2 checkpoint: `msa`, `lm300m` or `lm600m`. The language-model selections pass `--use-esm` automatically; confidence output depends on the upstream checkpoint and runtime version.                                                    |
| `--colabfold2_kernel_backend`    | `auto`  | Fused-kernel implementation used by `colabfold_batch`; `auto` is the upstream default, but on GPUs missing from the bundled JAX device table (e.g. L40S, V100) it fails with `No supported GPU devices found` and `cuda_legacy` must be set. |
| `--colabfold_num_recycles`       | `3`     | The number of times model outputs are recycled as iterative refinement input. This parameter is shared with the ColabFold mode.                                                                                                              |

> You can override any of these parameters via the command line or a params file.
