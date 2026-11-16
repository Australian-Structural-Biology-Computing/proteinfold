#!/usr/bin/env python3
"""Unit tests for the MSA adaptor scripts.

Run with:  python3 tests/test_msa_adaptors.py

These exercise bin/af2_msa_to_colabfold_json.py (AlphaFold2 MSA pipeline
output -> AlphaFold3-format input JSON, including the species-pairing logic
ported from alphafold's pipeline_multimer) and bin/check_af3_msa_json.py
(pass-through validation) against synthetic fixtures that mimic jackhmmer /
hhblits output. No external dependencies; Python >= 3.6.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF2_SCRIPT = os.path.join(PROJECT, "bin", "af2_msa_to_colabfold_json.py")
CHECK_SCRIPT = os.path.join(PROJECT, "bin", "check_af3_msa_json.py")
JSON2FASTA_SCRIPT = os.path.join(PROJECT, "bin", "alphafold3_json_to_af2_fasta.py")


def write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(content)
    return path


class Af2AdaptorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.msas = os.path.join(self.tmp, "msas")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    # --- fixtures -------------------------------------------------------

    def _heteromer_fixture(self):
        """Chain A (10 aa) and chain B (12 aa) sharing one uniprot species."""
        chain_map = {
            "A": {"sequence": "MKTAYIAKQR", "description": "chain_A"},
            "B": {"sequence": "MKTAYIAKQRAA", "description": "chain_B"},
        }
        write(
            os.path.join(self.msas, "chain_id_map.json"),
            json.dumps(chain_map, indent=4),
        )

        # Query row carries a gap column; other rows' residues there must
        # become lowercase insertions in the a3m output.
        write(
            os.path.join(self.msas, "A", "uniref90_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "A0/1-10              MKTAY-IAKQR\n"
            "tr|Q1AAAB|Q1AAAB_HUMAN/1-10 MKTAYiIAKQR\n"
            "tr|Q3CCCC|Q3CCCC_ECOLI/1-9  MKTAYiIAKQ-\n"
            "//\n",
        )
        write(
            os.path.join(self.msas, "A", "bfd_uniref_hits.a3m"),
            ">101 BFD query\nMKTAY-IAKQR\n"
            ">tr|Q4DDDD|Q4DDDD_MOUSE/1-10\nMKTAY-IADQR\n"
            ">tr|Q1AAAB|Q1AAAB_HUMAN/1-10\nMKTAYiIAKQR\n",  # duplicate
        )
        write(
            os.path.join(self.msas, "A", "mgnify_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "A0/1-10              MKTAY-IAKQR\n"
            "MGNNNNNN/1-8         MKTAYiIAK-\n"
            "//\n",
        )
        write(
            os.path.join(self.msas, "A", "uniprot_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "A0/1-10              MKTAY-IAKQR\n"
            "tr|Q1AAAB|Q1AAAB_HUMAN/1-100  MKTAYiIAKQR\n"
            "tr|Q3CCCC|Q3CCCC_ECOLI/1-90   MKTAYiIAK-\n"
            "//\n",
        )

        write(
            os.path.join(self.msas, "B", "uniref90_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "B0/1-12              MKTAYIAKQRAA\n"
            "tr|Q5EEEE|Q5EEEE_HUMAN/1-12 MKTAYIAKQRAA\n"
            "tr|Q3CCCC|Q3CCCC_ECOLI/1-10 MKTAYIAKQR--\n"
            "//\n",
        )
        write(
            os.path.join(self.msas, "B", "bfd_uniref_hits.a3m"),
            ">101 BFD query\nMKTAYIAKQRAA\n",
        )
        write(
            os.path.join(self.msas, "B", "mgnify_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "B0/1-12              MKTAYIAKQRAA\n"
            "MGNNNNNN/1-10        MKTAYIAKQ--\n"
            "//\n",
        )
        write(
            os.path.join(self.msas, "B", "uniprot_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "B0/1-12              MKTAYIAKQRAA\n"
            "tr|Q5EEEE|Q5EEEE_HUMAN/1-100  MKTAYIAKQRAA\n"
            "tr|Q3CCCC|Q3CCCC_ECOLI/1-110  MKTAYIAKQR--\n"
            "//\n",
        )

        fasta = write(
            os.path.join(self.tmp, "sample.fasta"),
            ">SAMPLE chain one\nMKTAYIAKQR\n"
            "> SAMPLE chain two\nMKTAYIAKQRAA\n",
        )
        return fasta

    def _run(self, fasta, msa_dir, output, name, extra=()):
        return subprocess.run(
            [sys.executable, AF2_SCRIPT, "--fasta", fasta, "--msa_dir", msa_dir,
             "--output", output, "--name", name] + list(extra),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )

    # --- tests ----------------------------------------------------------

    def test_heteromer_pairing_and_unpaired_blocks(self):
        fasta = self._heteromer_fixture()
        out = os.path.join(self.tmp, "out.json")
        result = self._run(fasta, self.msas, out, "SAMPLE")
        self.assertEqual(result.returncode, 0, result.stderr.decode())

        job = json.load(open(out))
        self.assertEqual(job["dialect"], "alphafold3")
        self.assertEqual(job["version"], 2)
        self.assertEqual(job["modelSeeds"], [1])

        prot_a, prot_b = [s["protein"] for s in job["sequences"]]
        self.assertEqual(prot_a["id"], ["A"])
        self.assertEqual(prot_b["id"], ["B"])
        self.assertEqual(prot_a["sequence"], "MKTAYIAKQR")

        # unpaired: query first, insertions lowercase, bfd/mgnify rows
        # appended, rows aligned-identical to the query deduplicated
        # (matching alphafold's make_msa_features dedup).
        unpaired = prot_a["unpairedMsa"]
        self.assertTrue(unpaired.startswith(">A0/1-10\nMKTAYIAKQR\n"))
        self.assertNotIn("Q1AAAB_HUMAN", unpaired)     # aligned-duplicate of query
        self.assertIn("Q3CCCC_ECOLI", unpaired)        # uniref90 hit
        self.assertIn("Q4DDDD_MOUSE", unpaired)        # bfd hit
        self.assertIn("MGNNNNNN", unpaired)            # mgnify hit

        # paired: equal depth, aligned synthetic descriptions, row 0 = query
        rows_a = [l for l in prot_a["pairedMsa"].splitlines() if l.startswith(">")]
        rows_b = [l for l in prot_b["pairedMsa"].splitlines() if l.startswith(">")]
        self.assertEqual(rows_a, rows_b)
        self.assertEqual(rows_a[0], ">101")
        self.assertEqual(rows_a[1], ">tr|AAAAAA|AAAAAA_00001")
        # HUMAN row paired, then ECOLI row (product-of-indices ordering)
        self.assertIn("MKTAYiIAKQR", prot_a["pairedMsa"])
        self.assertIn("MKTAYIAKQRAA", prot_b["pairedMsa"])
        self.assertIn("MKTAYiIAK-", prot_a["pairedMsa"])
        self.assertIn("MKTAYIAKQR--", prot_b["pairedMsa"])

    def test_duplicate_chain_shares_msa(self):
        fasta = self._heteromer_fixture()
        # Add chain C identical to A: no dedicated msa dir exists for it.
        chain_map = json.load(open(os.path.join(self.msas, "chain_id_map.json")))
        chain_map["C"] = {"sequence": "MKTAYIAKQR", "description": "chain_C"}
        write(
            os.path.join(self.msas, "chain_id_map.json"),
            json.dumps(chain_map, indent=4),
        )
        fasta = write(
            os.path.join(self.tmp, "triple.fasta"),
            ">c1\nMKTAYIAKQR\n>c2\nMKTAYIAKQRAA\n>c3\nMKTAYIAKQR\n",
        )
        out = os.path.join(self.tmp, "triple.json")
        result = self._run(fasta, self.msas, out, "TRIPLE")
        self.assertEqual(result.returncode, 0, result.stderr.decode())

        prot_a, prot_b, prot_c = [s["protein"] for s in json.load(open(out))["sequences"]]
        self.assertEqual([p["id"] for p in (prot_a, prot_b, prot_c)],
                         [["A"], ["B"], ["C"]])
        self.assertEqual(prot_a["unpairedMsa"], prot_c["unpairedMsa"])
        self.assertEqual(prot_a["pairedMsa"], prot_c["pairedMsa"])

    def test_monomer_has_no_paired_msa(self):
        fasta = self._heteromer_fixture()
        # Keep only chain A.
        chain_map = json.load(open(os.path.join(self.msas, "chain_id_map.json")))
        write(
            os.path.join(self.msas, "chain_id_map.json"),
            json.dumps({"A": chain_map["A"]}, indent=4),
        )
        fasta = write(os.path.join(self.tmp, "mono.fasta"), ">mono\nMKTAYIAKQR\n")
        out = os.path.join(self.tmp, "mono.json")
        result = self._run(fasta, self.msas, out, "MONO")
        self.assertEqual(result.returncode, 0, result.stderr.decode())

        protein = json.load(open(out))["sequences"][0]["protein"]
        self.assertEqual(protein["pairedMsa"], "")
        self.assertTrue(protein["unpairedMsa"])

    def test_homomer_reuses_unpaired_as_paired(self):
        fasta = self._heteromer_fixture()
        chain_map = json.load(open(os.path.join(self.msas, "chain_id_map.json")))
        write(
            os.path.join(self.msas, "chain_id_map.json"),
            json.dumps({"A": chain_map["A"],
                        "B": {"sequence": "MKTAYIAKQR", "description": "homo_B"}},
                       indent=4),
        )
        fasta = write(os.path.join(self.tmp, "homo.fasta"),
                      ">h1\nMKTAYIAKQR\n>h2\nMKTAYIAKQR\n")
        out = os.path.join(self.tmp, "homo.json")
        result = self._run(fasta, self.msas, out, "HOMO")
        self.assertEqual(result.returncode, 0, result.stderr.decode())

        prot_a, prot_b = [s["protein"] for s in json.load(open(out))["sequences"]]
        self.assertTrue(prot_a["pairedMsa"])
        self.assertEqual(prot_a["pairedMsa"], prot_b["pairedMsa"])

    def test_missing_uniprot_for_multimer_fails_loudly(self):
        fasta = self._heteromer_fixture()
        os.remove(os.path.join(self.msas, "B", "uniprot_hits.sto"))
        out = os.path.join(self.tmp, "bad.json")
        result = self._run(fasta, self.msas, out, "BAD")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"uniprot_hits.sto", result.stderr)

    def test_aligned_and_species_helpers(self):
        # Sanity check of the helpers used by dedup/pairing.
        from importlib.util import spec_from_file_location, module_from_spec
        spec = spec_from_file_location("af2_adaptor", AF2_SCRIPT)
        module = spec_from_file_location and module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.aligned_sequence("MKTAY-IAKQR"), "MKTAY-IAKQR")
        self.assertEqual(module.aligned_sequence("MKTAYiIAKQR"), "MKTAYIAKQR")
        self.assertEqual(module.species_identifier("tr|Q1AAAB|Q1AAAB_HUMAN/1-10"), "HUMAN")
        self.assertEqual(module.species_identifier("A0/1-10"), "")


class CheckAf3MsaJsonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _job(self, **chain_fields):
        protein = {"id": ["A"], "sequence": "MKTAYIAKQR",
                   "modifications": [], "templates": []}
        protein.update(chain_fields)
        return {"dialect": "alphafold3", "version": 2, "name": "test",
                "sequences": [{"protein": protein}], "modelSeeds": [1]}

    def _run(self, job):
        src = write(os.path.join(self.tmp, "in.json"), json.dumps(job))
        dst = os.path.join(self.tmp, "out.json")
        result = subprocess.run(
            [sys.executable, CHECK_SCRIPT, src, "--output", dst],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        return result, dst

    def test_valid_job_passes_through(self):
        result, dst = self._run(self._job(unpairedMsa=">101\nMKTAYIAKQR\n",
                                          pairedMsa=""))
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertTrue(os.path.exists(dst))

    def test_missing_msa_fields_rejected(self):
        result, dst = self._run(self._job())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"neither unpairedMsa nor pairedMsa", result.stderr)
        self.assertFalse(os.path.exists(dst))

    def test_paired_only_rejected(self):
        result, dst = self._run(self._job(pairedMsa=">101\nMKTAYIAKQR\n"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"no unpairedMsa", result.stderr)

    def test_multi_job_document_rejected(self):
        job = self._job(unpairedMsa=">101\nMKTAYIAKQR\n")
        result, dst = self._run([job, job])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"exactly one", result.stderr)


class JsonToFastaTests(unittest.TestCase):
    """alphafold3_json_to_af2_fasta.py: protein chains -> multi-record FASTA."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, job_or_list, name="SAMPLE"):
        src = write(os.path.join(self.tmp, "in.json"), json.dumps(job_or_list))
        dst = os.path.join(self.tmp, "out.fasta")
        result = subprocess.run(
            [sys.executable, JSON2FASTA_SCRIPT, src, "--id", name, "--output", dst],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        return result, dst

    def test_expands_copies_and_preserves_order(self):
        job = {
            "dialect": "alphafold3", "version": 2, "name": "x",
            "sequences": [
                {"protein": {"id": ["A"], "sequence": "MKTAYIAKQR"}},
                {"rna": {"id": ["R1"], "sequence": "GGCC"}},
                {"protein": {"id": ["B", "C"], "sequence": "GGGHIJKLMN"}},
                {"ligand": {"ccdCodes": ["ATP"]}},
            ],
            "modelSeeds": [1],
        }
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        with open(dst) as handle:
            content = handle.read()
        records = [line[1:] for line in content.splitlines() if line.startswith(">")]
        sequences = [
            line for line in content.splitlines()
            if line and not line.startswith(">")
        ]
        self.assertEqual(records, ["SAMPLE_1", "SAMPLE_2", "SAMPLE_3"])
        self.assertEqual(sequences, ["MKTAYIAKQR", "GGGHIJKLMN", "GGGHIJKLMN"])

    def test_rejects_protein_free_json(self):
        job = {
            "dialect": "alphafold3", "version": 2, "name": "x",
            "sequences": [{"ligand": {"ccdCodes": ["ATP"]}}],
            "modelSeeds": [1],
        }
        result, dst = self._run(job)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"at least one protein entity", result.stderr)

    def test_rejects_multi_job_document(self):
        job = {"sequences": [{"protein": {"id": ["A"], "sequence": "MKT"}}]}
        result, dst = self._run([job, job])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"one AlphaFold3 job per samplesheet row", result.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
