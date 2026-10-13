import pathlib
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
from axion_fopt_next import Catalogue, AxionPotential
from axion_fopt_next.interchange import write_arrival_map


@unittest.skipUnless(shutil.which("c++"), "C++ compiler required for adapter check")
class InterchangeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.probe = pathlib.Path(cls.tmp.name) / "probe"
        root = pathlib.Path(__file__).resolve().parents[1]
        subprocess.run(["c++", "-std=c++20", "-Wall", "-Wextra", "-Werror",
                        str(root / "adapters/cosmolattice/probe.cpp"),
                        "-o", str(cls.probe)], check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def probe_map(self, path, stage, origin=.5, omega=2.):
        return subprocess.run([str(self.probe), str(path), str(stage), str(origin),
                               str(omega), ".4", ".2", "1.5"], capture_output=True, text=True)

    def test_one_bubble_units_and_nonmonotonic_stages(self):
        arrival = Catalogue(4., .7, [1.], [[.2, 1., 2.]]).lattice_arrival_times(4)
        path = pathlib.Path(self.tmp.name) / "bubble.afopt"
        write_arrival_map(path, arrival, 4.)
        potential = AxionPotential(.2, 1.5)
        # Include reversal in stage time; the background must have no time history.
        for stage in (0., 8., 2., 6., 1.):
            result = self.probe_map(path, stage)
            self.assertEqual(result.returncode, 0, result.stderr)
            actual = np.loadtxt(result.stdout.splitlines())
            tau = .5 + stage / 2.
            expected = np.column_stack([
                potential.mass_squared(tau, arrival).ravel() / 4.,
                potential.force(.4, tau, arrival).ravel() / 4.,
                potential.energy(.4, tau, arrival).ravel() / 4.])
            np.testing.assert_allclose(actual, expected, rtol=2e-15)

    def test_boundary_and_never_converted(self):
        arrival = np.array([1., np.inf, 2., 3., 4., 5., 6., 7.]).reshape(2, 2, 2)
        path = pathlib.Path(self.tmp.name) / "boundary.afopt"
        write_arrival_map(path, arrival, 2.)
        result = self.probe_map(path, 1.)  # tau=1 exactly
        self.assertEqual(result.returncode, 0, result.stderr)
        values = np.loadtxt(result.stdout.splitlines())
        np.testing.assert_allclose(values[:, 0], [1.5**2/4] + [.2**2/4]*7)

    def test_malformed_maps_and_units(self):
        path = pathlib.Path(self.tmp.name) / "malformed.afopt"
        write_arrival_map(path, np.ones((2, 2, 2)), 2.)
        good = path.read_bytes()
        for bad in (b"wrong", good[:-1], good + b"x",
                    good[:8] + (2**32-1).to_bytes(4, "little") + good[12:],
                    good[:20] + np.array([np.nan], dtype="<f8").tobytes() + good[28:]):
            path.write_bytes(bad)
            self.assertNotEqual(self.probe_map(path, 0.).returncode, 0)
        path.write_bytes(good)
        self.assertNotEqual(self.probe_map(path, 0., omega=0.).returncode, 0)
        self.assertNotEqual(self.probe_map(path, np.nan).returncode, 0)


class WriterTests(unittest.TestCase):
    def test_invalid_maps(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "bad.afopt"
            for array in (np.ones((2, 3, 2)), np.full((1, 1, 1), np.nan),
                          np.full((1, 1, 1), -np.inf)):
                with self.assertRaises(ValueError): write_arrival_map(path, array, 1.)
            self.assertFalse(path.exists())
