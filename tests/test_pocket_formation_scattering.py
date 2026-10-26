"""Continuum scattering checks against exact constant/step solutions."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from formation_scattering_reference import scatter,step_amplitudes


class FormationScattering(unittest.TestCase):
    def test_constant_mass_is_transparent(self):
        row=scatter(.6,kappa=0.,extent=8.)
        self.assertLess(row['reflected_flux'],1e-20)
        self.assertAlmostEqual(row['transmitted_flux'],1.,delta=2e-10)

    def test_sharp_step_matches_amplitude_and_flux(self):
        for omega in (.5,1.25):
            row=scatter(omega,profile='step',extent=8.,rtol=2e-12,atol=2e-14)
            amplitude,r,t=step_amplitudes(omega,.1,1.)
            np.testing.assert_allclose(complex(row['reflection_real'],row['reflection_imag']),amplitude,rtol=2e-8,atol=2e-9)
            self.assertAlmostEqual(row['reflected_flux'],r,delta=2e-9)
            self.assertAlmostEqual(row['transmitted_flux'],t,delta=2e-9)

    def test_subthreshold_and_input_rejection(self):
        row=scatter(.6782081956648878)
        self.assertAlmostEqual(row['reflected_flux'],1.,delta=1e-10)
        self.assertEqual(row['transmitted_flux'],0.)
        self.assertLess(row['flux_residual'],1e-10)
        with self.assertRaises(ValueError):scatter(.05)
