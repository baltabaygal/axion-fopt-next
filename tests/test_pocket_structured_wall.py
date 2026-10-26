"""Independent checks of the proposed structured potential reference."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
import structured_wall as wall


class StructuredWall(unittest.TestCase):
    def test_force_and_curvature_are_potential_derivatives(self):
        phi=np.linspace(-.03,1.03,60);eps=1e-6
        gradient=(wall.potential(phi+eps,2.,.01,12.)-wall.potential(phi-eps,2.,.01,12.))/(2*eps)
        curvature=(wall.force(phi+eps,2.,.01,12.)-wall.force(phi-eps,2.,.01,12.))/(2*eps)
        np.testing.assert_allclose(gradient,wall.force(phi,2.,.01,12.),atol=2e-7,rtol=2e-7)
        np.testing.assert_allclose(curvature,wall.curvature(phi,2.,.01,12.),atol=2e-7,rtol=2e-7)
        self.assertEqual(wall.potential(0.,2.,.01,12.),.01)
        self.assertEqual(wall.potential(1.,2.,.01,12.),0.)
        self.assertEqual(wall.force(0.,2.,.01,12.),0.)
        self.assertEqual(wall.force(1.,2.,.01,12.),0.)
        self.assertAlmostEqual(wall.curvature(0.,2.,.01,12.),13.)
        self.assertAlmostEqual(wall.curvature(1.,2.,.01,12.),13.12)

    def test_default_kink_and_tension_match_analytic_solution(self):
        z=np.linspace(-8.,8.,101);profile=wall.kink(2.,0.)
        np.testing.assert_allclose(profile(z),.5*(1+np.tanh(z/2)),atol=3e-13)
        self.assertAlmostEqual(wall.tension(2.,0.),1/6,places=13)

    def test_chosen_tilted_potential_has_only_one_barrier(self):
        x=np.linspace(1e-6,1-1e-6,10001);f=wall.force(x,2.,.01,12.)
        self.assertEqual(np.count_nonzero((f[:-1]>0)&(f[1:]<0)),1)
        self.assertEqual(np.count_nonzero((f[:-1]<0)&(f[1:]>0)),0)
        for invalid in (-1.,16.,np.nan):
            with self.assertRaises(ValueError):
                wall.potential(.3,2.,.01,invalid)


if __name__=='__main__':
    unittest.main()
