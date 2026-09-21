# ProteinFold test fixtures

Small, deterministic inputs and CPU-stage outputs used by the ASBC nf-test
suite live here. Model weights and databases are deliberately not committed.

Real CPU tests and GPU tests backed entirely by the shared databases use the
existing `s3://proteinfold-dataset/test-data/mini_dbs` CI dataset by default.
Set `PROTEINFOLD_TEST_DB` to the root of a compatible local microDB mirror to
override the S3 dataset.

GPU input fixtures derived from CPU tests are stored under `cpu_outputs/` so
the inference tests run independently. Their provenance, generating database
identity and SHA-256 checksums are recorded below. They are refreshed only
when the corresponding CPU output contract intentionally changes.

## GPU test inputs

The inference tests use the checked-in CPU outputs by default. The defaults
can be overridden for local validation:

- `PROTEINFOLD_AF2_FEATURES`: validated `RUN_ALPHAFOLD2_MSA` `features.pkl`.
- `PROTEINFOLD_AF3_DATA_JSON`: validated `RUN_ALPHAFOLD3_DATAPIPELINE`
  `*_data.json`.
- `PROTEINFOLD_AF3_PARAMS`: licensed AlphaFold3 `af3.bin`; defaults to the
  protected test-data bucket and is never committed.

The currently accepted CPU artifacts were generated from the canonical nf-core
test-datasets `T1026.fasta`
with ProteinFold commit `3e884c97b512a897c8294a1e1d788dd050206fc4`, random
seed `0`, and the shared ProteinFold microDB dataset:

| CPU artifact                    | SHA-256                                                            |
| ------------------------------- | ------------------------------------------------------------------ |
| AlphaFold2 `T1026_features.pkl` | `51bd8ba32a548c5a76e2f1fa0ce5ab0f38b7383baa0f5999dda6f2f1de5b513b` |
| AlphaFold3 `T1026_data.json`    | `b987a743a617ae18d8ef337efdb15cfda9cc953acbd57d0b6319caaffb5cafe1` |

Each GPU test verifies this checksum before submitting its inference process.
When the corresponding CPU contract changes, validate the new artifact first,
then update the documented checksum and GPU test together.
