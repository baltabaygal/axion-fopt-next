import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from planar_initial import prepare,PlanarSeed
from planar_boundary_initial import extended_seed,moving_final_core


class PlanarBoundary(unittest.TestCase):
    def test_common_core_initial_data_are_exact(self):
        old_x,old=prepare(4096)
        x,state,raw=extended_seed(.125)
        core=(old_x>=32)&(old_x<=288)
        overlap=(x>=32)&(x<=288)
        np.testing.assert_array_equal(x[overlap],old_x[core])
        np.testing.assert_array_equal(state[overlap],old[core])
        self.assertTrue(np.all(state[:20,[1,3]]==0.))
        self.assertTrue(np.all(state[-20:,[1,3]]==0.))
        self.assertLess(np.max(abs(state[:,[1,3]]-raw)),1e-10)

    def test_moving_wall_core_is_preserved(self):
        old_x,old=prepare(4096,PlanarSeed(speed=.3))
        x,state,_=extended_seed(.125,speed=.3)
        np.testing.assert_array_equal(state[(x>=32)&(x<=288)],old[(old_x>=32)&(old_x<=288)])

    def test_final_comparison_excludes_swept_second_wall(self):
        lo,hi=moving_final_core(.3,180.)
        self.assertGreater(lo,.3*180.)
        self.assertLess(hi,256-.3*180.)
        self.assertLess(lo,146.)
        self.assertGreater(hi,146.)
        with self.assertRaises(ValueError):moving_final_core(.3,500.)

    def test_strong_packet_core_is_preserved(self):
        old_x,old=prepare(4096,PlanarSeed(speed=.3,amplitude=.16))
        x,state,_=extended_seed(.125,speed=.3,amplitude=.16)
        np.testing.assert_array_equal(state[(x>=32)&(x<=288)],old[(old_x>=32)&(old_x<=288)])
