import sys,unittest
from pathlib import Path
from dataclasses import replace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from planar_distributed_initial import prepare,audit,SlabSeed

class DistributedSeed(unittest.TestCase):
    def test_symmetric_walls_and_stress(self):
        row=audit(4096)
        self.assertLess(row['reflection_symmetry_error'],1e-13)
        self.assertLess(row['axion_antisymmetry_error'],1e-13)
        self.assertEqual(row['tv_energy_per_area'],0.)
        self.assertAlmostEqual(row['modes']['linear_phase_average_pressure'],.0002,places=15)
    def test_empty_matches_walls(self):
        _,s,_=prepare(2048)
        _,empty,_=prepare(2048,replace(SlabSeed(),pressure_ratio=0))
        np.testing.assert_array_equal(s[:,[0,2]],empty[:,[0,2]])
        self.assertTrue(np.all(empty[:,[1,3]]==0))
    def test_same_continuum_modes_across_grids(self):
        rows=[audit(n) for n in (2048,4096,8192)]
        self.assertTrue(all(r['modes']==rows[0]['modes'] for r in rows))
        e=[r['axion_energy_per_area'] for r in rows]
        self.assertGreater(abs(e[0]-e[1])/abs(e[1]-e[2]),12)
