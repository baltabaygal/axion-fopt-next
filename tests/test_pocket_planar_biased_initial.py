import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from planar_biased_initial import audit,prepare_biased
from planar_initial import prepare,PlanarSeed
from paper_reference import Parameters

class BiasedSeed(unittest.TestCase):
    def test_force_sign_and_normalization(self):
        for dv in (.0001,.0002,.001):
            row=audit(2048,p=Parameters(dv=dv))
            self.assertLess(row['relative_force_error'],1e-8)
            self.assertLess(row['thin_wall_initial_acceleration'],0)
    def test_profile_is_explicitly_approximate(self):
        seed=PlanarSeed(speed=.3,amplitude=.16)
        with self.assertRaises(ValueError):prepare(2048,seed,Parameters(dv=.0002))
        _,s=prepare_biased(2048,seed,Parameters(dv=.0002))
        _,reference=prepare(2048,seed,Parameters(dv=0))
        np.testing.assert_array_equal(s,reference)
    def test_axion_energy_includes_gradients(self):
        row=audit()
        self.assertAlmostEqual(row['initial_axion_energy_per_area'],.091458590884567,places=10)
