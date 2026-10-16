"""Necessary resolution bounds and coordinate-unit independence."""
import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('pq_feasibility',Path(__file__).resolve().parents[1]/'checks/pq_strings/feasibility.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

class PQFeasibilityTests(unittest.TestCase):
    def test_target_bound_and_relaxed_reference(self):
        self.assertEqual(audit.bounds(3.)['minimum_N'],36450)
        self.assertEqual(audit.bounds(3.,core_cells=1.,mass_ratio=.2)['minimum_N'],7290)
        for n in (32,64,128):self.assertLess(audit.best_ratio(n),1.1)

    def test_current_setup_and_unit_independence(self):
        case=audit.evaluate(64,24.,20.,27.,.08,8.,3.)
        self.assertAlmostEqual(case['wall_cells'],.2962962962962963)
        c=.5
        scaled=audit.evaluate(64,24*c,20*c,27*c,.08/c**2,8.,3.)
        for key in ('ma_ms','wall_cells','core_cells','kmin_tau_star','L_over_tau_end'):
            self.assertAlmostEqual(case[key],scaled[key])

    def test_bound_expresses_mass_and_core_interval(self):
        b=audit.bounds(3.);r=3.;B=b['minimum_L_over_tau_star'];sigma_min=r**4/.16
        self.assertAlmostEqual(4*B*r*sigma_min,b['core_bound_N'])
        sigma_max=128/(4*B*r)
        self.assertGreater(sigma_min,sigma_max)
