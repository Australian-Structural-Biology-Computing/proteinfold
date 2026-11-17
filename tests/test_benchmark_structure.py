"""Run with the OpenStructure container: python3 -m unittest discover -s tests -p test_benchmark_structure.py."""

import json
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "benchmark_structure.py"


def structure(chains):
    lines = []
    for chain in chains:
        for residue in range(1, 13):
            x = residue * 3.8
            for name, dx, y, z, element in [
                ("N", 0.0, 0.0, 0.0, "N"),
                ("CA", 1.3, 0.7, 0.0, "C"),
                ("C", 2.5, 0.0, 0.0, "C"),
                ("O", 2.8, -1.2, 0.0, "O"),
                ("CB", 1.3, 0.7, 1.5, "C"),
            ]:
                lines.append(
                    f"ATOM  {len(lines) + 1:5d} {name:4s} ALA {chain}{residue:4d}    "
                    f"{x + dx:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}  \n"
                )
    return "".join(lines) + "END\n"


class BenchmarkStructureTest(unittest.TestCase):
    def run_comparison(self, model_chains, reference_chains):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            model = root / "model.pdb"
            reference = root / "reference.pdb"
            output = root / "scores.json"
            model.write_text(structure(model_chains))
            reference.write_text(structure(reference_chains))
            process = subprocess.run(
                ["python3", str(SCRIPT), "--model", str(model), "--reference", str(reference), "--output", str(output)],
                capture_output=True,
                text=True,
            )
            result = json.loads(output.read_text()) if output.exists() else None
            return process, result

    def test_rejects_multiple_chains(self):
        process, result = self.run_comparison("AB", "A")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("exactly one protein chain", process.stderr)
        self.assertIsNone(result)

    def test_score_validation(self):
        spec = importlib.util.spec_from_file_location("benchmark_structure", SCRIPT)
        benchmark = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(benchmark)
        data = {
            "status": "SUCCESS", "model_chains": ["A"], "reference_chains": ["B"],
            "chain_mapping": {"B": "A"}, "usalign_mapping": {"B": "A"}, "lddt": 1.0, "tm_score": 0.8,
        }
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "scores.json"

            def fake_run(command, check):
                Path(command[command.index("-o") + 1]).write_text(json.dumps(data))

            args = ["benchmark_structure.py", "--model", "model.pdb", "--reference", "ref.pdb", "--output", str(output)]
            with patch.object(sys, "argv", args), patch.object(benchmark, "check_protein"), patch.object(
                benchmark.subprocess, "run", side_effect=fake_run
            ):
                benchmark.main()
                self.assertEqual(json.loads(output.read_text())["lddt"], 1.0)
                data["tm_score"] = None
                with self.assertRaisesRegex(ValueError, "invalid tm_score"):
                    benchmark.main()


if __name__ == "__main__":
    unittest.main()
