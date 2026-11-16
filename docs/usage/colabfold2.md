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

## Compiled-kernel cache

JIT compilation dominates the wall clock of short predictions: every task that
sees a new input shape pays the full compile cost again. `--colabfold2_compile_cache`
points the pipeline at a persistent directory that survives between tasks and
runs, so a shape compiles once per machine:

```bash
nextflow run nf-core/proteinfold \
    --colabfold2_compile_cache /scratch/share/colabfold_compile_cache \
    ...
```

Two mechanisms share the directory:

- The AlphaFold2 backends keep colabfold's own per-shape executable cache
  (`--compile-cache`). Entries are keyed by GPU, driver and package versions
  and by the model config and input shapes, and are only written when fused
  kernels are enabled (`--use-fast-kernels`, e.g. via
  `--colabfold2_kernel_backend pallas`); the bundled `colabfold_compile`
  utility can pre-fill or prune the directory (`colabfold_compile --clean`
  removes entries left behind by upgrades).
- The AlphaFold3-family backends compile once per sequence-length bucket; the
  directory is exported as `JAX_COMPILATION_CACHE_DIR` so JAX's persistent
  compilation cache reuses those compiled kernels across tasks and runs.

The directory is bound into the inference containers automatically (the
pipeline creates it on the host at launch if it does not exist yet), so no
custom `--singularity-args` are needed. It is deliberately _not_ staged into
the task directories: staged inputs are hashed, and a cache that grows over
time would then invalidate `-resume` on every run. This way tasks hash only
the directory's path, and the cache contents can evolve freely between
resumed runs.

Both mechanisms cache compiled executables that a run loads and executes, so
a cache directory is code, not data: anyone who can write to it can run code
with the privileges of every user whose runs consume it. Keep the cache
per-user (e.g. `~/.cache/proteinfold/compile`, mode `0700`) and only share it
within a group of mutually trusting accounts that runs under the same
account boundary.
The directory must be writable inside the container. If your container engine
does not mount host paths by default, bind it explicitly, e.g.
`--singularity-args '--bind /scratch/share/colabfold_compile_cache'`.
## Batching

Both expensive stages can process several samples in a single task:

| Parameter                | Default | Description                                                                                          |
| ------------------------ | ------- | ---------------------------------------------------------------------------------------------------- |
| `--msa_batch_size`       | `-1`    | Number of samples processed together during MSA search; `-1` runs all inputs in one batch.           |
| `--inference_batch_size` | `1`     | Number of samples processed together during structure prediction; `-1` runs all inputs in one batch. |

Batching structure prediction groups several samples into one `colabfold_batch` task so the model weights are loaded (and, where applicable, compiled) once for the whole batch instead of once per sample. Batched inference requires one query per samplesheet row (use `--split_fasta` for multi-record FASTA inputs) and is not supported with `--use_msa_server`; batched outputs are published flat under the selected mode token (files keep their per-sample names) instead of in per-sample subdirectories.

## File Structure

If `--colabfold2_params_path` is not set and no weights cache is available under a local `--db`/`--colabfold_db`, the pipeline downloads the required weights in a separate Nextflow process. The weights root has this structure:

<details markdown="1">
<summary>Directory structure</summary>

```
<colabfold2_params_path>/params/
├── alphafold_params*.npz                               (alphafold2 / alphafold2_ptm)
├── alphafold_multimer_params*.npz                      (alphafold2_multimer_v3)
├── download_finished.txt                               (alphafold2_ptm marker)
├── download_complexes_multimer_v3_finished.txt         (alphafold2_multimer_v3 marker)
└── af3/
    ├── libcifpp/
    │   └── components.cif
    ├── esmc_300m|esmc_600m/
    │   └── <variant>.bin.zst
    └── <model>[-int8]/
        ├── <model>.bin.zst
        ├── <model>.lm.npz
        └── download_<model>_<precision>_finished.txt
```

The `alphafold2` modes use ColabFold's own AlphaFold2 parameter layout: `<root>/params/` with the `alphafold_params_*.npz` files and success markers written by `colabfold.download.download_alphafold_params` (the same markers its downloader skips on). When only `--db`/`--colabfold_db` is given, the implicit `<colabfold_db>/colabfold2_models` cache is checked for those markers too, for the AlphaFold2 parameter sets exactly as for the AlphaFold3-family weights — a complete cache avoids the download entirely; a partial one triggers a fresh task-local download of everything.

Downloaded weights are consumed directly from the task work directory and are not published into the run outputs.

Chai1 and Boltz2 use revision-qualified completion markers instead: `download_chai1_<precision>_r3_finished.txt` and `download_boltz2_<precision>_r1_finished.txt`.

</details>

`<model>` is the selected upstream model, suffixed with `-int8` when `--colabfold2_weights_precision int8` is set; the `<model>.lm.npz` checkpoint and matching `esmc_300m` or `esmc_600m` trunk are only required by the ESMFold2 language-model selections.

## Additional Arguments

See the [ColabFold2](https://github.com/sokrypton/ColabFold) documentation for a full description of additional arguments. The arguments supported by the proteinfold workflow are described briefly below:

| Parameter                        | Default | Description                                                                                                                                                                                                                                                                                                                                                                                                      |
| -------------------------------- | ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--colabfold2_weights_precision` | `fp32`  | Checkpoint precision: `fp32` or `int8`. Official AlphaFold3 weights support `fp32` only, and the AlphaFold2 models provide no int8 weights.                                                                                                                                                                                                                                                                      |
| `--esmfold2_model`               | `msa`   | ESMFold2 checkpoint: `msa`, `lm300m` or `lm600m`. The language-model selections pass `--use-esm` automatically; confidence output depends on the upstream checkpoint and runtime version.                                                                                                                                                                                                                        |
| `--colabfold2_kernel_backend`    | `auto`  | Fused-kernel implementation used by `colabfold_batch`; `auto` is the upstream default, but on GPUs missing from the bundled JAX device table (e.g. L40S, V100) it fails with `No supported GPU devices found` and `cuda_legacy` must be set.                                                                                                                                                                     |
| `--colabfold_num_recycles`       | `3`     | The number of times model outputs are recycled as iterative refinement input. This parameter is shared with the ColabFold mode.                                                                                                                                                                                                                                                                                  |
| `--colabfold2_compile_cache`     | none    | Persistent compiled-kernel cache directory, shared across tasks and runs so JIT compilation is paid once per machine. Passed to `colabfold_batch` as `--compile-cache` (effective for the AlphaFold2 backends together with `--use-fast-kernels`, see `--colabfold2_kernel_backend`) and exported as `JAX_COMPILATION_CACHE_DIR`, which the AlphaFold3-family backends use for their per-bucket JIT compilation. |

> You can override any of these parameters via the command line or a params file.
