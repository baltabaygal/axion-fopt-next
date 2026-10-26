"""Pocket-model potential, initial-state and I/O contracts, plus optional compiled-backend validation."""
import importlib.util
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
try:
    import scipy
except ImportError:
    scipy = None
if scipy is not None:
    spec = importlib.util.spec_from_file_location('pocket_reference', ROOT/'checks/pockets/model.py')
    pk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pk)


@unittest.skipIf(scipy is None, 'Pocket references require optional scipy')
class PocketContracts(unittest.TestCase):
    def test_force_is_potential_gradient(self):
        rng = np.random.default_rng(5)
        phi, a = rng.uniform(-.5, 1.5, (2, 100));eps = 1e-6
        numerical = ((pk.potential(phi+eps, a, 2., .01, 4.)-pk.potential(phi-eps, a, 2., .01, 4.))/(2*eps),
                     (pk.potential(phi, a+eps, 2., .01, 4.)-pk.potential(phi, a-eps, 2., .01, 4.))/(2*eps))
        for n, f in zip(numerical, pk.force(phi, a, 2., .01, 4.)):
            np.testing.assert_allclose(f, n, atol=2e-8, rtol=2e-8)

    def test_vacua_and_axion_mass(self):
        self.assertEqual(pk.dilaton_potential(1., 2., .01), 0.)
        self.assertAlmostEqual(pk.dilaton_potential(0., 2., .01), .01, places=15)
        self.assertEqual(pk.force(0., 0., 2., .01, 4.), (0., 0.))
        self.assertEqual(pk.force(1., 0., 2., .01, 4.), (0., 0.))
        # Axion is massless in false vacuum and has mass^2 = m_TV^2 in true vacuum.
        self.assertEqual(pk.force(0., .3, 2., .01, 4.)[1], 0.)
        self.assertAlmostEqual(pk.force(1., .3, 2., .01, 4.)[1], 1.2)

    def test_quartic_coupling_and_false_vacuum_curvature(self):
        rng=np.random.default_rng(11);phi,a=rng.uniform(-.5,1.5,(2,60));eps=1e-6
        numerical=((pk.potential(phi+eps,a,2.,.01,4.,4)-pk.potential(phi-eps,a,2.,.01,4.,4))/(2*eps),
                   (pk.potential(phi,a+eps,2.,.01,4.,4)-pk.potential(phi,a-eps,2.,.01,4.,4))/(2*eps))
        for observed,expected in zip(pk.force(phi,a,2.,.01,4.,4),numerical):
            np.testing.assert_allclose(observed,expected,rtol=2e-8,atol=2e-8)
        for amplitude in (0.,1.,10.):
            curvature=(pk.force(eps,amplitude,2.,.01,4.,4)[0]-pk.force(-eps,amplitude,2.,.01,4.,4)[0])/(2*eps)
            self.assertAlmostEqual(curvature,1.,delta=1e-8)
        self.assertAlmostEqual(pk.force(1.,.3,2.,.01,4.,4)[1],1.2)
        for invalid in (0,2.5,3,6):
            with self.assertRaises(ValueError):
                pk.force(.3,.2,2.,.01,4.,invalid)

    def test_gas_field_is_grid_independent(self):
        modes = pk.gas_modes(16., .5, 1.3, 2)
        coarse, _ = pk.gas_field(modes, 12);fine, _ = pk.gas_field(modes, 24)
        k = 2*np.pi*modes['vectors']/16.
        for field, n in ((coarse, 12), (fine, 24)):
            x = pk.grid(n, 16.)
            point = np.array([x[1], x[3], x[5]])
            direct = (modes['field']*np.exp(1j*(k@point))).sum().real
            self.assertAlmostEqual(field[1, 3, 5], direct, places=10)
        with self.assertRaises(ValueError):
            pk.gas_field(modes, 4)

    def test_gas_amplitude_hits_target_energy(self):
        modes = pk.gas_modes(16., .4, 1., 1)
        amp = pk.gas_amplitude(3., 16, 16., 5., 1., modes, 4., 1., 4., 8., .01)
        s = pk.pocket_state(16, 16., 5., 1., modes, amp, 4., 1.)
        e = pk.energies(s, 16., 8., .01, 4.)
        self.assertAlmostEqual(e['axion_kinetic']+e['axion_gradient']+e['coupling'], 3., places=10)
        self.assertGreater(e['axion_fv'], .9*3.)

    def test_state_roundtrip_and_validation(self):
        s = pk.pocket_state(6, 6., 2., .5)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'state.bin'
            pk.write_initial(path, s, 6.)
            # PKIN0001 header: magic, n, box. Snapshots (PKST0001) add the time u.
            raw = path.read_bytes();self.assertEqual(raw[:8], b'PKIN0001');self.assertEqual(len(raw), 20+32*6**3)
            path.write_bytes(b'PKST0001'+struct.pack('<I2d', 6, 6., .5)+raw[20:])
            meta, back = pk.read_state(path)
            self.assertEqual((meta['n'], meta['box'], meta['u']), (6, 6., .5))
            self.assertTrue(np.array_equal(back, s))
            path.write_bytes(path.read_bytes()[:-8])
            with self.assertRaises(ValueError):
                pk.read_state(path)
            with self.assertRaises(ValueError):
                pk.write_initial(path, s[..., :3], 6.)

    def test_radial_modes_reach_reflective_sphere_limit(self):
        r0 = 5.
        modes = pk.radial_modes(r0, 1., 400., .3, 2.5, 6., h=.01, profile=lambda r: (r > r0).astype(float))
        low = [m for m in modes if m['l'] <= 1 and m['n'] <= 1]
        self.assertEqual(len(low), 4)
        for m in low:
            self.assertAlmostEqual(m['omega']*r0/pk.bessel_zeros(m['l'], 2)[m['n']], 1., delta=3e-3)
            self.assertAlmostEqual(np.trapezoid(m['u']**2, m['r']), 1., places=6)

    def test_eigen_gas_is_confined_and_grid_independent(self):
        modes = pk.radial_modes(4., 1., 3., .5, 1.2, 7.)
        self.assertGreater(len(modes), 0)
        coarse, coarse_p = pk.eigen_gas(modes, 16, 16., 4)
        fine, fine_p = pk.eigen_gas(modes, 48, 16., 4)
        # Refinement by 3 keeps every coarse cell centre: fine index 3i+1.
        np.testing.assert_allclose(fine[1::3, 1::3, 1::3], coarse, atol=1e-12)
        np.testing.assert_allclose(fine_p[1::3, 1::3, 1::3], coarse_p, atol=1e-12)
        s = pk.eigen_state(48, 16., 4., 1., modes, 1., 4)
        r = pk.radius(48, 16.);dx = 16/48
        grad = sum((np.roll(s[..., 1], -1, ax)-s[..., 1])**2 for ax in range(3))/(2*dx*dx)
        density = .5*s[..., 3]**2+grad+.5*9*s[..., 0]**2*s[..., 1]**2
        self.assertLess(density[r > 7].sum()/density.sum(), 1e-4)

    def test_gas_chunking_preserves_the_same_realisation_bitwise(self):
        modes=pk.radial_modes(4.,1.,3.,.5,1.2,7.)
        small=pk.eigen_gas(modes,24,16.,4,7.,chunk_size=13)
        full=pk.eigen_gas(modes,24,16.,4,7.,chunk_size=24**3)
        for left,right in zip(small,full):
            self.assertTrue(np.array_equal(left,right))

    def test_thin_wall_brackets(self):
        sigma = pk.tension(pk.lambda_for_width(2.))
        b = pk.brackets(20., sigma, .01, 3*(.01+2*sigma/17)*.85**4, 2756.7)
        self.assertAlmostEqual(b['adiabatic'], 17., places=8)
        self.assertGreater(b['energy_conserving'], 15.)


EXECUTABLE = os.environ.get('AFOPT_POCKET_EXECUTABLE')


@unittest.skipUnless(EXECUTABLE and scipy is not None, 'Set AFOPT_POCKET_EXECUTABLE to run backend validation')
class PocketBackend(unittest.TestCase):
    def test_validation_suite(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run([sys.executable, str(ROOT/'checks/pockets/validate.py'), '--output',
                                     str(Path(tmp)/'validation'), '--executable', EXECUTABLE],
                                    cwd=ROOT/'checks/pockets', capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout[-3000:]+result.stderr[-3000:])


if __name__ == '__main__':
    unittest.main()
