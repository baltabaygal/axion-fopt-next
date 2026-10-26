"""Cell-centre alignment in the periodic grid-refinement comparison."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from refine_planar_scattering import coarse_sample


class PlanarRefinement(unittest.TestCase):
    def test_fourier_shift_matches_common_physical_coordinates(self):
        box=8.;nc=32;nf=128
        def field(n):
            x=(np.arange(n)+.5)*box/n;k=2*np.pi/box
            return np.column_stack((np.cos(k*x),np.sin(3*k*x),np.cos(5*k*x+.3),np.ones(n)))
        np.testing.assert_allclose(coarse_sample(field(nf),nc,box),field(nc),rtol=0.,atol=5e-15)
        with self.assertRaises(ValueError):coarse_sample(field(nf),30,box)
