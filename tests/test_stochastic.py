"""Independent survival references and identical-event exclusion checks."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
from axion_fopt_next import sample_candidates,exclude_converted_births,sample_catalogue,Catalogue


class CandidateTests(unittest.TestCase):
    def test_same_events_and_phantom_equivalence(self):
        for seed in (0,1,42):
            candidates=sample_candidates([0.,1.],[4.,4.],3.,.5,seed)
            accepted=exclude_converted_births(candidates)
            direct=sample_catalogue([0.,1.],[4.,4.],3.,.5,seed)
            np.testing.assert_array_equal(accepted.times,direct.times)
            np.testing.assert_array_equal(accepted.positions,direct.positions)
            probes=np.random.default_rng(910).uniform(0,3,(512,3))
            np.testing.assert_allclose(candidates.arrival_times(probes),accepted.arrival_times(probes),atol=1e-12,rtol=0)
            for i,(birth,centre) in enumerate(zip(accepted.times,accepted.positions)):
                earlier=Catalogue(3.,.5,accepted.times[:i],accepted.positions[:i])
                self.assertGreater(earlier.arrival_times(centre),birth)


@unittest.skipUnless(importlib.util.find_spec("scipy"),"SciPy required for independent analytic reference")
class AnalyticTests(unittest.TestCase):
    def test_static_closed_form_and_cosmic_change_of_variables(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"checks/stochastic"))
        from analytic import CASES,survival,cosmic_survival,tabulated_survival
        for case in CASES:
            times=np.linspace(0,case.end,21)
            np.testing.assert_allclose(survival(case,times),cosmic_survival(case,times),atol=1e-9,rtol=0)
        static=CASES[0]; times=np.linspace(0,1,21)
        np.testing.assert_allclose(tabulated_survival([0,1],[4,4],times,.5),survival(static,times),atol=1e-13,rtol=0)
