"""Analytic hierarchy and fixed-a backreaction checks, without lattice runs."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from formation_hierarchy import candidate,false_branch
from paper_reference import Parameters


class FormationHierarchy(unittest.TestCase):
    def test_candidate_boost_and_grid_estimates(self):
        coarse=candidate();fine=candidate(n=256)
        self.assertTrue(all(coarse['gates'].values()))
        self.assertAlmostEqual(coarse['mass_true'],1.)
        self.assertFalse(coarse['repeated_head_on_estimates'][1]['encounters'][1]['below_exterior_threshold'])
        self.assertEqual(coarse['dirichlet_angular_mode_count'],8)
        self.assertEqual([v['ell'] for v in coarse['dirichlet_band_modes']],[1,2])
        gamma=1/np.sqrt(1-.3**2)
        for rest,lab in zip(coarse['wall_frame_penetration_lengths'],coarse['boosted_lab_penetration_cells']):
            self.assertAlmostEqual(lab,rest/gamma/.25)
        self.assertAlmostEqual(fine['initial_tanh_10_90_cells']/coarse['initial_tanh_10_90_cells'],256/192)
        self.assertGreater(coarse['long_lived_target_time'],coarse['initial_interface_to_image_gap'])

    def test_large_boost_is_not_automatically_reflective(self):
        row=candidate(speed=.999)
        self.assertFalse(row['gates']['initial_head_on_reflection'])
        self.assertFalse(row['gates']['boosted_lab_penetration_four_cells'])
        with self.assertRaises(ValueError):candidate(speed=1.)

    def test_shift_and_effective_curvature(self):
        p=Parameters(dv=.05);zero=false_branch(0.,p)
        self.assertEqual(zero['phi_minimum'],0.)
        for a in (2.,12.):
            row=false_branch(a,p);self.assertLess(row['phi_minimum'],0.)
            self.assertGreater(row['phi_curvature'],0.)
            self.assertLess(abs(row['phi_force_residual']),1e-12)
            h=1e-3
            numerical=(false_branch(a+h,p)['potential']-2*row['potential']+false_branch(a-h,p)['potential'])/h**2
            self.assertAlmostEqual(numerical,row['adiabatic_axion_curvature'],delta=2e-8)
        self.assertLess(false_branch(12.,p)['full_hessian_min_eigenvalue'],0.)
