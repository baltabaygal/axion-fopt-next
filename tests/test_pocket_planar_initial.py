"""Directionality and boost signs for coupled planar initial fields."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from planar_initial import PlanarSeed,prepare,audit,modes
from paper_reference import Parameters


class PlanarInitial(unittest.TestCase):
    def test_packet_right_moving_and_localized(self):
        row=audit(2048)
        self.assertLess(row['continuum_left_moving_fraction'],1e-25)
        self.assertLess(row['lattice_left_moving_fraction'],1e-6)
        self.assertLess(row['initial_energy_in_tv_fraction'],1e-10)
        self.assertLess(row['initial_energy_near_walls_fraction'],1e-10)

    def test_same_continuum_coefficients_on_refined_grids(self):
        seed=PlanarSeed();p=Parameters(dv=0.)
        _,coarse=prepare(2048,seed,p);x,fine=prepare(4096,seed,p)
        i,k,w,c=modes(seed,p)
        expected=sum(seed.amplitude*weight*np.cos(kk*(x-seed.centre)) for kk,weight in zip(k,c))
        np.testing.assert_array_equal(fine[:,1],expected)
        self.assertEqual(coarse.shape,(2048,4))
        self.assertLess(audit(4096)['continuum_kink_discrete_residual'],audit(2048)['continuum_kink_discrete_residual']/3.)

    def test_contracting_false_region_has_positive_wall_momentum(self):
        _,s=prepare(2048,PlanarSeed(speed=.3))
        self.assertGreater(s[0,2],0.)
        self.assertGreater(s[1024,2],0.)
        np.testing.assert_allclose(np.roll(s[:,0],1024),1-s[:,0],rtol=0.,atol=1e-13)
        np.testing.assert_allclose(np.roll(s[:,2],1024),s[:,2],rtol=0.,atol=1e-13)
        with self.assertRaises(ValueError):prepare(2048,p=Parameters(dv=.01))
