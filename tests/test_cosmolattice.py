"""Optional integration suite; enable with AFOPT_COSMOLATTICE_EXECUTABLE."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from axion_fopt_next.interchange import write_arrival_map

EXECUTABLE = os.environ.get("AFOPT_COSMOLATTICE_EXECUTABLE")


@unittest.skipUnless(EXECUTABLE, "Build CosmoLattice adapter and set AFOPT_COSMOLATTICE_EXECUTABLE")
class CosmoLatticeTests(unittest.TestCase):
    def test_field_evolution_references(self):
        script = Path(__file__).resolve().parents[1] / "examples/validate_cosmolattice.py"
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "validation"
            result = subprocess.run([sys.executable, str(script), "--executable", EXECUTABLE,
                                     "--output", str(output)], capture_output=True, text=True, timeout=180)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            summary = json.loads((output / "validation.json").read_text())
            self.assertEqual(len(summary["checks"]), 8)

    def test_reject_unsupported_or_invalid_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "empty.afopt"
            write_arrival_map(path, np.full((4,4,4), np.inf), 4.)
            defaults = {"N":4, "lSide":4, "dt":.01, "tMax":.1,
                        "expansion":"false", "evolver":"PV2", "arrival_map":path,
                        "result_dir": root / "result"}
            for override in ({"evolver":"LF"}, {"expansion":"true", "fixedBackground":"false"},
                             {"lSide":3}, {"omega_star":0}, {"dt":1},
                             {"dt":.03}, {"theta_initial":"nan"}, {"mass_true":-1}, {"theta_initial":"garbage"}, {"N":4.5}):
                params = defaults | override
                result = subprocess.run([EXECUTABLE, *[f"{k}={v}" for k,v in params.items()]],
                                        capture_output=True, text=True, timeout=30)
                self.assertNotEqual(result.returncode, 0, override)
                self.assertFalse((root / "result").exists(), override)
