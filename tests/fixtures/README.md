# ProteinFold test fixtures

Small, deterministic inputs and CPU-stage outputs used by the ASBC nf-test
suite live here. Model weights and databases are deliberately not committed.

Real CPU and GPU tests use the database root selected by
`PROTEINFOLD_TEST_DB`. On Katana this should be:

```text
/srv/scratch/sbf-pipelines/databases/proteinfold_microdbs
```

GPU input fixtures derived from CPU tests must include a provenance file with
the generating ProteinFold commit, command, database identity, parameters and
SHA-256 checksum. They are refreshed only when the corresponding CPU output
contract intentionally changes.

## GPU test inputs

The inference tests require explicit paths so they cannot silently select an
old workflow result:

- `PROTEINFOLD_AF2_FEATURES`: validated `RUN_ALPHAFOLD2_MSA` `features.pkl`.
- `PROTEINFOLD_AF3_DATA_JSON`: validated `RUN_ALPHAFOLD3_DATAPIPELINE`
  `*_data.json`.
- `PROTEINFOLD_AF3_PARAMS`: licensed AlphaFold3 `af3.bin`; never committed.

The currently accepted CPU artifacts were generated from `ubiquitin.fasta`
with ProteinFold commit `3e884c97b512a897c8294a1e1d788dd050206fc4`, random
seed `0`, and the Katana microDB root shown above:

| CPU artifact                     | SHA-256                                                            |
| -------------------------------- | ------------------------------------------------------------------ |
| AlphaFold2 `features.pkl`        | `212cb501978f2a93fab1c91b071a8329d92cc763efeb69b08a47caf6ad2449e5` |
| AlphaFold3 `ubiquitin_data.json` | `1e4503b48ffd630e5d6f99054bd5bc729267bb4c62ad39b341e2e4707f6fa562` |

Each GPU test verifies this checksum before submitting its inference process.
When the corresponding CPU contract changes, validate the new artifact first,
then update the documented checksum and GPU test together.
