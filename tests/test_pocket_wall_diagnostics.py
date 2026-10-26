"""Independent geometrical checks for local wall-width diagnostics."""
import sys
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from wall_diagnostics import local_widths,sampled_gradient


class WallDiagnostics(unittest.TestCase):
    @staticmethod
    def legacy_gradient(phi,points,dx):
        gradients=[(np.roll(phi,-1,a)-np.roll(phi,1,a))/(2*dx) for a in range(phi.ndim)]
        return np.stack([map_coordinates(g,points.T,order=1,mode='grid-wrap') for g in gradients],-1)

    def test_sampled_gradients_match_original_across_periodic_cuts(self):
        rng=np.random.default_rng(87)
        for ndim in (2,3):
            phi=rng.normal(size=(12,)*ndim)
            points=rng.uniform(-3,15,size=(400,ndim))
            np.testing.assert_allclose(sampled_gradient(phi,points,.2),
                                       self.legacy_gradient(phi,points,.2),atol=1e-13,rtol=0.)

    def test_unchanged_profile_measurements(self):
        for ndim in (2,3):
            x=(np.arange(24)+.5)*.5-6
            coords=np.meshgrid(*([x]*ndim),indexing='ij')
            r=np.sqrt(sum(c*c for c in coords))
            for phi in (.5*(1+np.tanh((r-3)/.7)),.2+.3*(1+np.tanh((r-3)/.7))):
                observed=local_widths(phi,.5,reach=2.)
                with patch('wall_diagnostics.sampled_gradient',side_effect=self.legacy_gradient):
                    expected=local_widths(phi,.5,reach=2.)
                for key in expected:
                    if key in ('min_cells','median_cells') and expected[key] is not None:
                        self.assertAlmostEqual(observed[key],expected[key],delta=1e-11)
                    else:
                        self.assertEqual(observed[key],expected[key])

    def test_spherical_tanh_has_analytic_normal_thickness(self):
        # A curved wall verifies that widths follow normals rather than a
        # coordinate axis (whose width diverges near a tangent).
        n, dx, radius, width = 128, .25, 10., 1.
        x=(np.arange(n)+.5)*dx-n*dx/2
        r=np.hypot(x[:,None],x[None,:])
        phi=.5*(1+np.tanh((r-radius)/width))
        result=local_widths(phi,dx,reach=4.)
        expected=2*np.arctanh(.8)*width/dx
        self.assertEqual(result['ambiguous'],0)
        self.assertGreater(result['contour_edges'],100)
        self.assertAlmostEqual(result['min_cells'],expected,delta=.2)
        self.assertAlmostEqual(result['median_cells'],expected,delta=.2)

    def test_missing_thresholds_fail_without_dropping_edges(self):
        x=(np.arange(64)+.5)*.5-16
        r=np.hypot(x[:,None],x[None,:])
        phi=.2+.6*.5*(1+np.tanh((r-8)))
        result=local_widths(phi,.5)
        self.assertGreater(result['contour_edges'],0)
        self.assertEqual(result['ambiguous'],result['contour_edges'])
        self.assertEqual(result['measured'],0)
        self.assertEqual(result['ambiguity_counts'],dict(missing_threshold_bracket=result['ambiguous']))
        self.assertIsNone(result['min_cells'])
        self.assertFalse(result['diagnostic_pass'])

    def test_no_pocket_does_not_pass_resolution(self):
        result=local_widths(np.ones((16,16)),.5)
        self.assertEqual(result['contour_edges'],0)
        self.assertFalse(result['diagnostic_pass'])


if __name__=='__main__':
    unittest.main()
