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
JSON2A3M_SCRIPT = os.path.join(PROJECT, "bin", "af3_msa_json_to_a3m.py")


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

    def test_monomer_layout_without_chain_map(self):
        """The monomer data pipeline writes alignments at the msa dir root."""
        write(
            os.path.join(self.msas, "uniref90_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "A0/1-10              MKTAY-IAKQR\n"
            "tr|Q3CCCC|Q3CCCC_ECOLI/1-9  MKTAYiIAKQ-\n"
            "//\n",
        )
        write(
            os.path.join(self.msas, "bfd_uniref_hits.a3m"),
            ">101 BFD query\nMKTAY-IAKQR\n",
        )
        write(
            os.path.join(self.msas, "mgnify_hits.sto"),
            "# STOCKHOLM 1.0\n\n"
            "A0/1-10              MKTAY-IAKQR\n"
            "MGNNNNNN/1-8         MKTAYiIAK-\n"
            "//\n",
        )
        fasta = write(
            os.path.join(self.tmp, "sample.fasta"),
            ">SAMPLE one chain\nMKTAYIAKQR\n",
        )
        out = os.path.join(self.tmp, "out.json")
        result = self._run(fasta, self.msas, out, "SAMPLE")
        self.assertEqual(result.returncode, 0, result.stderr.decode())

        job = json.load(open(out))
        prot = job["sequences"][0]["protein"]
        self.assertEqual(prot["id"], ["A"])
        self.assertEqual(prot["sequence"], "MKTAYIAKQR")
        self.assertTrue(prot["unpairedMsa"].startswith(">A0/1-10\nMKTAYIAKQR\n"))
        self.assertIn("Q3CCCC_ECOLI", prot["unpairedMsa"])
        self.assertEqual(prot["pairedMsa"], "")

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


class JsonToA3mTests(unittest.TestCase):
    """af3_msa_json_to_a3m.py: AF3-format MSA JSON -> ColabFold a3m."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _job(self, unpaired_a, paired_a, unpaired_b=None, paired_b=None):
        sequences = [{"protein": {"id": ["A"], "sequence": "MKTAY",
                                  "unpairedMsa": unpaired_a, "pairedMsa": paired_a}}]
        if unpaired_b is not None:
            sequences.append({"protein": {"id": ["B"], "sequence": "GGGC",
                                          "unpairedMsa": unpaired_b, "pairedMsa": paired_b}})
        return {"dialect": "alphafold3", "version": 2, "name": "t",
                "sequences": sequences, "modelSeeds": [1]}

    def _run(self, job):
        src = write(os.path.join(self.tmp, "in.json"), json.dumps(job))
        dst = os.path.join(self.tmp, "out.a3m")
        result = subprocess.run(
            [sys.executable, JSON2A3M_SCRIPT, src, "--output", dst],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        return result, dst

    def test_complex_with_paired_msa(self):
        job = self._job(
            ">101\nMKTAY\n>102\nMKTAiY\n>103\nMKTAY\n",
            ">101\nMKTAY\n>tr|x|x_00001\nMKTAY\n",
            ">101\nGGGC\n>104\nGGGC\n",
            ">101\nGGGC\n>tr|x|x_00001\nGG-G-\n",
        )
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        lines = open(dst).read().splitlines()
        # header: chain lengths and cardinalities
        self.assertEqual(lines[0], "#5,4\t1,1")
        # query row: tab-joined descriptions, concatenated sequences
        self.assertEqual(lines[1], ">101\t101")
        self.assertEqual(lines[2], "MKTAYGGGC")
        # paired rows aligned by row index, tab-joined descriptions
        self.assertEqual(lines[3], ">tr|x|x_00001\ttr|x|x_00001")
        self.assertEqual(lines[4], "MKTAYGG-G-")
        # unpaired block: per-chain rows padded to full complex width
        self.assertEqual(lines[6], "MKTAY----")
        self.assertEqual(lines[12], "-----GGGC")

    def test_monomer_unpaired_only(self):
        job = self._job(">101\nMKTAY\n>102\nMKTAiY\n", "")
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        lines = open(dst).read().splitlines()
        self.assertEqual(lines[0], "#5\t1")
        self.assertEqual(lines[1], ">101")
        self.assertEqual(lines[2], "MKTAY")
        self.assertIn("MKTAiY", lines)

    def test_crlf_msa_rows_are_normalized(self):
        job = self._job(">101\r\nMKTAY\r\n>102\r\nMKTAiY\r\n", "")
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        with open(dst, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0], "#5\t1")
        self.assertEqual(lines[1:], [">101", "MKTAY", ">102", "MKTAiY"])

    def test_homodimer_msa_width_matches_cardinality(self):
        job = self._job(">101\nMKTAY\n>102\nMKTAiY\n", "")
        job["sequences"][0]["protein"]["id"] = ["A", "B"]
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        with open(dst, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0], "#5\t2")
        self.assertEqual([len(line.replace("i", "")) for line in lines[1:] if not line.startswith(">")], [10, 10, 10, 10])

    def test_separate_identical_entities_share_one_msa_block(self):
        job = self._job(
            ">101\nMKTAY\n>102\nMKTAiY\n",
            ">101\nMKTAY\n>tr|x|x_00001\nMKTAY\n",
        )
        job["sequences"].append({"protein": {
            "id": ["B"], "sequence": "MKTAY",
            "unpairedMsa": ">101\nMKTAY\n>ignored\nAAAAA\n",
            "pairedMsa": ">101\nMKTAY\n>ignored\nAAAAA\n",
        }})
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        with open(dst, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0], "#5\t2")
        self.assertNotIn("ignored", "\n".join(lines))
        self.assertEqual({len(line.replace("i", "")) for line in lines[1:] if not line.startswith(">")}, {10})

    def test_unequal_paired_depths_rejected(self):
        job = self._job(
            ">101\nMKTAY\n>102\nMKTAY\n", ">101\nMKTAY\n>tr|a|a_00001\nMKTAY\n",
            ">101\nGGGC\n", ">101\nGGGC\n>tr|a|a_00001\nGGGCC\n>tr|b|b_00002\nGGGCC\n",
        )
        result, dst = self._run(job)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"unequal depth", result.stderr)

    def test_missing_unpaired_reports_absent_msa(self):
        job = self._job("", "")
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 2)
        self.assertIn(b"no unpairedMsa", result.stderr)

    def test_protein_free_json_reports_absent_msa(self):
        # A ligand-only JSON has non-protein entities, which are reported
        # distinctly (exit 3) before the missing-protein check applies.
        job = {"dialect": "alphafold3", "version": 2, "name": "t",
               "sequences": [{"ligand": {"ccdCodes": ["ATP"]}}], "modelSeeds": [1]}
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 3)
        self.assertIn(b"non-protein entities", result.stderr)

    def test_non_protein_entities_rejected_distinctly(self):
        job = {"dialect": "alphafold3", "version": 2, "name": "t",
               "sequences": [
                   {"protein": {"id": ["A"], "sequence": "MKTAY",
                                 "unpairedMsa": ">101\nMKTAY\n"}},
                   {"ligand": {"id": ["L1"], "ccdCodes": ["ATP"]}},
               ], "modelSeeds": [1]}
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 3)
        self.assertIn(b"non-protein entities", result.stderr)
        self.assertIn(b"ligand", result.stderr)

    def test_malformed_protein_remains_an_error(self):
        job = {"dialect": "alphafold3", "version": 2, "name": "t",
               "sequences": [{"protein": {"id": ["A"]}}], "modelSeeds": [1]}
        result, dst = self._run(job)
        self.assertEqual(result.returncode, 1)
        self.assertIn(b"missing its sequence", result.stderr)


STUB_TREE = """alphafold/__init__.py:
    pass
alphafold/data/__init__.py:
    pass
alphafold/data/parsers.py:
    class TemplateHit:
        def __init__(self, name, aligned_cols, hit_sequence, indices_query,
                     indices_hit, query):
            self.name = name
            self.aligned_cols = aligned_cols
            self.hit_sequence = hit_sequence
            self.indices_query = indices_query
            self.indices_hit = indices_hit
            self.query = query

    def convert_stockholm_to_a3m(sto, remove_first_row_gaps=False):
        return sto

    def parse_hhr(hhr_string):
        hits = []
        records = [r for r in hhr_string.split(">") if r.strip()]
        for index, record in enumerate(records, start=1):
            lines = record.splitlines()
            name, sequence = lines[0], "".join(lines[1:])
            name = name.split("/")[0]
            ungapped = sequence.replace("-", "")
            hits.append(TemplateHit(
                name=name,
                aligned_cols=len(ungapped),
                hit_sequence=sequence,
                indices_query=list(range(len(ungapped))),
                indices_hit=list(range(len(ungapped))),
                query="",
            ))
        return hits

    def parse_hmmsearch_a3m(query_sequence, a3m_string, skip_first=False):
        hits = []
        records = [r for r in a3m_string.split(">") if r.strip()]
        for index, record in enumerate(records, start=1):
            lines = record.lstrip(">").splitlines()
            name, sequence = lines[0], "".join(lines[1:])
            # like alphafold's _parse_hmmsearch_description: pdb0_A/1-4 -> pdb0_A
            name = name.split("/")[0]
            ungapped = sequence.replace("-", "")
            hits.append(TemplateHit(
                name=name,
                aligned_cols=len(ungapped),
                hit_sequence=sequence,
                indices_query=list(range(len(ungapped))),
                indices_hit=list(range(len(ungapped))),
                query=query_sequence,
            ))
        return hits
alphafold/data/mmcif_parsing.py:
    class _Result:
        def __init__(self, header):
            self.mmcif_object = type("Obj", (), {"header": header})()

    def parse(file_id, mmcif_string):
        # Test fixtures carry their release date on a "# release_date:" line.
        header = {"release_date": None}
        for line in mmcif_string.splitlines():
            if line.startswith("# release_date:"):
                header["release_date"] = line.split(":", 1)[1].strip()
        return _Result(header)
colabfold/__init__.py:
    pass
colabfold/alphafold3/__init__.py:
    pass
colabfold/alphafold3/templates.py:
    class _Template:
        def __init__(self, mmcif, mapping):
            self.mmcif = mmcif
            self.query_to_template_map = mapping

    def build_templates(hits, query_sequence, cif_dir, max_templates=20):
        return [_Template("cif-" + hit.name, {0: 0, 1: 1}) for hit in hits]
"""


def _make_stub_tree(root):
    """Materialise stub alphafold/colabfold packages for the subprocess."""
    import re

    for match in re.finditer(r"^(\S+\.py):\n((?:    .*\n|\n)*)", STUB_TREE, re.M):
        path, body = match.group(1), match.group(2)
        path = os.path.join(root, path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        body = "\n".join(line[4:] if line.startswith("    ") else line
                          for line in body.splitlines())
        with open(path, "w") as handle:
            handle.write(body + "\n")


class Af2TemplatesToJsonTests(unittest.TestCase):
    """Template conversion in bin/af2_msa_to_colabfold_json.py."""

    SCRIPT = AF2_SCRIPT

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        _make_stub_tree(os.path.join(self.tmp.name, "stubs"))
        with open(os.path.join(self.tmp.name, "obsolete.dat"), "w") as handle:
            handle.write("")

    def run_script(self, max_date="2100-01-01"):
        env = dict(os.environ)
        env["PYTHONPATH"] = os.path.join(self.tmp.name, "stubs")
        return subprocess.run(
            [
                sys.executable, self.SCRIPT,
                "--fasta", os.path.join(self.tmp.name, "sample.fasta"),
                "--msa_dir", os.path.join(self.tmp.name, "raw"),
                "--output", os.path.join(self.tmp.name, "out.json"),
                "--name", "T1026",
                "--template_mmcif_dir", os.path.join(self.tmp.name, "cifs"),
                "--obsolete", os.path.join(self.tmp.name, "obsolete.dat"),
                "--max-template-date", max_date,
            ],
            capture_output=True, text=True, env=env,
        )

    def _write_fixture(self, release_dates=("2020-01-01", "2020-01-01")):
        msa_dir = os.path.join(self.tmp.name, "raw", "T1026", "msas")
        os.makedirs(os.path.join(msa_dir, "A"))
        os.makedirs(os.path.join(msa_dir, "B"))
        chain_map = {
            "A": {"sequence": "MKVLAAAEEEKK", "description": "chainA"},
            "B": {"sequence": "GGWCCCFFFHHH", "description": "chainB"},
        }
        with open(os.path.join(msa_dir, "chain_id_map.json"), "w") as handle:
            json.dump(chain_map, handle)
        with open(os.path.join(msa_dir, "A", "pdb_hits.sto"), "w") as handle:
            handle.write(">pdb0_A/1-4\nMKVLAAAEEEKK\n")
        with open(os.path.join(msa_dir, "B", "pdb_hits.sto"), "w") as handle:
            handle.write(">pdb1_B/1-3\nGGACCCFFFHHH\n")
        for chain_id, sequence in (("A", "MKVLAAAEEEKK"), ("B", "GGWCCCFFFHHH")):
            with open(os.path.join(
                    msa_dir, chain_id, "uniref90_hits.sto"), "w") as handle:
                handle.write(
                    f"# STOCKHOLM 1.0\n\n{chain_id}0/1-12  {sequence}\n//\n")
            with open(os.path.join(
                    msa_dir, chain_id, "bfd_uniref_hits.a3m"), "w") as handle:
                handle.write(f">101 BFD query\n{sequence}\n")
            with open(os.path.join(
                    msa_dir, chain_id, "mgnify_hits.sto"), "w") as handle:
                handle.write(
                    f"# STOCKHOLM 1.0\n\n{chain_id}0/1-12  {sequence}\n//\n")
            with open(os.path.join(
                    msa_dir, chain_id, "uniprot_hits.sto"), "w") as handle:
                handle.write(
                    f"# STOCKHOLM 1.0\n\n{chain_id}0/1-12  {sequence}\n//\n")
        os.makedirs(os.path.join(self.tmp.name, "cifs"), exist_ok=True)
        for index, date in enumerate(release_dates):
            with open(os.path.join(
                    self.tmp.name, "cifs", f"pdb{index}.cif"), "w") as handle:
                handle.write(f"# release_date: {date}\ndata_x")
        with open(os.path.join(self.tmp.name, "sample.fasta"), "w") as handle:
            handle.write(">A\nMKVLAAAEEEKK\n>B\nGGWCCCFFFHHH\n")
        sample = {
            "dialect": "alphafold3", "version": 2, "name": "T1026",
            "sequences": [
                {"protein": {"id": ["A"], "sequence": "MKVLAAAEEEKK",
                             "templates": [], "unpairedMsa": ">101\nMKVLAAAEEEKK\n"}},
                {"protein": {"id": ["B"], "sequence": "GGWCCCFFFHHH",
                             "templates": [], "unpairedMsa": ">101\nGGWCCCFFFHHH\n"}},
            ],
        }
        with open(os.path.join(self.tmp.name, "sample.json"), "w") as handle:
            json.dump(sample, handle)

    def test_templates_populated_and_dated(self):
        self._write_fixture(release_dates=("2020-01-01", "2500-01-01"))
        with open(os.path.join(
                self.tmp.name, "raw", "T1026", "msas", "A",
                "pdb_hits.sto"), "w") as handle:
            # Second hit is a near-duplicate subsequence of the query.
            handle.write(">pdb0_A/1-4\nMKIAAAAEEEKK\n>pdb1_A/1-4\nMKVLAAAEEEKK\n")
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        document = json.load(open(os.path.join(self.tmp.name, "out.json")))
        templates_a = document["sequences"][0]["protein"]["templates"]
        templates_b = document["sequences"][1]["protein"]["templates"]
        # pdb1 is dropped as a duplicate; chain B's hit is after the cutoff.
        self.assertEqual(len(templates_a), 1)
        self.assertEqual(templates_a[0]["mmcif"], "cif-pdb0_A")
        self.assertEqual(templates_a[0]["queryIndices"], [0, 1])
        self.assertEqual(templates_a[0]["templateIndices"], [0, 1])
        self.assertEqual(templates_b, [])

    def test_obsolete_remap_and_removal(self):
        self._write_fixture()
        with open(os.path.join(
                self.tmp.name, "raw", "T1026", "msas", "A",
                "pdb_hits.sto"), "w") as handle:
            handle.write(
                ">pdbold_A/1-4\nMKIAAAAEEEKK\n>pdbgone_A/1-4\nMKIAAAAEEEKK\n>pdb0_A/1-4\nMKIEAAAEEEKK\n")
        with open(os.path.join(self.tmp.name, "obsolete.dat"), "w") as handle:
            handle.write(
                "OBSLTE 31-Jul-2020 pdb0_A pdbold_A\n"
                "OBSLTE 31-Jul-2020 pdbgone_A\n"
            )
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        document = json.load(open(os.path.join(self.tmp.name, "out.json")))
        templates = document["sequences"][0]["protein"]["templates"]
        # pdbold remaps onto pdb0 (deduplicated), pdbgone is fully removed.
        self.assertEqual(len(templates), 1)
        self.assertEqual(templates[0]["mmcif"], "cif-pdb0_A")

    def test_short_template_dropped(self):
        self._write_fixture()
        with open(os.path.join(
                self.tmp.name, "raw", "T1026", "msas", "A",
                "pdb_hits.sto"), "w") as handle:
            handle.write(">pdb0_A/1-4\nMKIAAAAEEEKK\n>pdb9_A/1-2\nMKIAA\n")
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        document = json.load(open(os.path.join(self.tmp.name, "out.json")))
        templates = document["sequences"][0]["protein"]["templates"]
        self.assertEqual(len(templates), 1)
        self.assertEqual(templates[0]["mmcif"], "cif-pdb0_A")

    def test_monomer_layout_with_hhr(self):
        """Monomer presets write pdb_hits.hhr at the msa dir root."""
        msa_dir = os.path.join(self.tmp.name, "raw", "T1026", "msas")
        os.makedirs(msa_dir)
        os.makedirs(os.path.join(self.tmp.name, "cifs"))
        with open(os.path.join(self.tmp.name, "cifs", "pdb0.cif"), "w") as handle:
            handle.write("# release_date: 2020-01-01\ndata_x")
        with open(os.path.join(msa_dir, "uniref90_hits.sto"), "w") as handle:
            handle.write("# STOCKHOLM 1.0\n\nA0/1-12  MKVLAAAEEEKK\n//\n")
        with open(os.path.join(msa_dir, "mgnify_hits.sto"), "w") as handle:
            handle.write("# STOCKHOLM 1.0\n\nA0/1-12  MKVLAAAEEEKK\n//\n")
        with open(os.path.join(msa_dir, "pdb_hits.hhr"), "w") as handle:
            handle.write(">pdb0_A/1-12\nMKIAAAAEEEKK\n")
        with open(os.path.join(self.tmp.name, "sample.fasta"), "w") as handle:
            handle.write(">A\nMKVLAAAEEEKK\n")
        sample = {
            "dialect": "alphafold3", "version": 2, "name": "T1026",
            "sequences": [
                {"protein": {"id": ["protA"], "sequence": "MKVLAAAEEEKK",
                             "templates": []}},
            ],
        }
        with open(os.path.join(self.tmp.name, "sample.json"), "w") as handle:
            json.dump(sample, handle)
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        document = json.load(open(os.path.join(self.tmp.name, "out.json")))
        templates = document["sequences"][0]["protein"]["templates"]
        self.assertEqual(len(templates), 1)
        self.assertEqual(templates[0]["mmcif"], "cif-pdb0_A")

    def test_missing_sto_leaves_empty(self):
        self._write_fixture()
        os.remove(os.path.join(
            self.tmp.name, "raw", "T1026", "msas", "B", "pdb_hits.sto"))
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        document = json.load(open(os.path.join(self.tmp.name, "out.json")))
        self.assertEqual(document["sequences"][1]["protein"]["templates"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
