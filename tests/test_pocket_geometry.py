"""Periodic translation checks for pocket position and stored energy."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from model import dilaton_sphere,radius
from moving_snapshot_budget import moving_budget
from pocket_geometry import periodic_centroid


class PocketGeometry(unittest.TestCase):
    def test_centroid_unwraps_a_translated_pocket(self):
        mask=radius(16,8.)<1.5
        shifted=np.roll(mask,(7,-6,5),axis=(0,1,2))
        centre,_=periodic_centroid(shifted,8.)
        np.testing.assert_allclose(centre,[7.5,1.,6.5],atol=1e-13)

    def test_uniform_or_absent_false_vacuum_has_no_defined_centroid(self):
        for mask in (np.ones((8,8,8),bool),np.zeros((8,8,8),bool)):
            with self.assertRaises(ValueError):
                periodic_centroid(mask,8.)

    def test_moving_energy_is_invariant_under_periodic_translation(self):
        n=16;box=16.;r=radius(n,box)
        state=np.zeros((n,n,n,4));state[...,0]=dilaton_sphere(n,box,4.,1.)
        x=(np.arange(n)+.5)-box/2
        state[...,1]=np.exp(-r*r/8)*np.cos(.7*x[:,None,None])
        state[...,3]=np.exp(-r*r/8)*np.sin(.7*x[:,None,None])
        shifted=np.roll(state,(5,-4,3),axis=(0,1,2))
        a=moving_budget(state,box,6.,8.,.01,4.,4)
        b=moving_budget(shifted,box,6.,8.,.01,4.,4)
        for key in ('stored_total_energy','stored_axion_energy','stored_phi_energy','global_total_energy'):
            self.assertAlmostEqual(a[key],b[key],delta=1e-12)
        np.testing.assert_allclose(a['stored_momentum'],b['stored_momentum'],atol=1e-12)


if __name__=='__main__':
    unittest.main()
