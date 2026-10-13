"""Independent checks of optional smooth-wall/quadratic adapter modes."""
import os
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from axion_fopt_next import Catalogue
from axion_fopt_next.interchange import write_arrival_map

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"checks/single_bubble"))
from reference import integrate, mass, observables, read_state

EXECUTABLE = os.environ.get("AFOPT_COSMOLATTICE_EXECUTABLE")


class FourierTests(unittest.TestCase):
    def test_zero_mode_and_plane_wave_normalization(self):
        # Uniform massive scalar with phi=1, pi=0 has n=1/2.
        stats, _ = observables(np.ones((8,8,8)),np.zeros((8,8,8)),2*np.pi)
        self.assertAlmostEqual(stats["number_continuum"],.5)
        self.assertAlmostEqual(stats["zero_mode_number"],.5)
        x = np.arange(8)*2*np.pi/8
        field = np.broadcast_to(np.cos(x)[:,None,None],(8,8,8))
        stats, _ = observables(field,np.zeros_like(field),2*np.pi)
        self.assertAlmostEqual(stats["number_continuum"],np.sqrt(2)/4)
        self.assertLess(stats["parseval_relative_error"],1e-14)


@unittest.skipUnless(EXECUTABLE,"Set AFOPT_COSMOLATTICE_EXECUTABLE")
class BubbleAdapterTests(unittest.TestCase):
    def test_irregular_arrival_cache_fallback(self):
        # No repeated arrivals: automatic cache fallback must preserve all cells.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            arrival=np.linspace(.1,2.,64).reshape(4,4,4)
            write_arrival_map(root/"map.afopt",arrival,4.)
            states=[]
            for cache in ("true","false"):
                out=root/cache
                args=dict(N=4,lSide=4,dt=.01,tMax=.5,expansion="false",evolver="PV2",
                          theta_initial=.7,mass_false=.2,mass_true=1.5,potential="cosine",
                          wall_width=.15,profile_wall_speed=.7,arrival_map=root/"map.afopt",
                          result_dir=out,snapshot_format="binary",mass_cache=cache,output_every=7)
                proc=subprocess.run([EXECUTABLE,*[f"{k}={v}" for k,v in args.items()]],capture_output=True,text=True,timeout=30)
                self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
                _,state=read_state(out/"final.bin")
                states.append(state)
            np.testing.assert_array_equal(*states)

    def test_smooth_quadratic_field_and_energy(self):
        if importlib.util.find_spec("scipy") is None:
            self.skipTest("SciPy required for independent PDE reference")
        from scipy.integrate import trapezoid
        n, box, v, end = 8, 2*np.pi, .5, 2.
        arrival = Catalogue(box,v,[0.],[[box/2]*3]).lattice_arrival_times(n)
        _, reference = integrate(arrival,box,v,end)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_arrival_map(root/"map.afopt",arrival,box)
            errors = []
            momentum_errors = []
            for dt in (.02,.01):
                out = root/f"dt{dt}"
                args = dict(N=n,lSide=box,dt=dt,tMax=end,expansion="false",evolver="PV2",
                            theta_initial=1.,theta_prime_initial=0.,tau_origin=0.,omega_star=1.,
                            mass_false=0.,mass_true=1.,potential="quadratic",wall_width=.15,
                            profile_wall_speed=v,arrival_map=root/"map.afopt",result_dir=out,
                            snapshot_format="binary",energy_diagnostics="true",output_every=1)
                proc = subprocess.run([EXECUTABLE,*[f"{k}={val}" for k,val in args.items()]],capture_output=True,text=True,timeout=60)
                self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
                meta, initial = read_state(out/"initial.bin")
                _, final = read_state(out/"final.bin")
                np.testing.assert_array_equal(initial[...,0],1.)
                np.testing.assert_array_equal(initial[...,1],0.)
                np.testing.assert_array_equal(initial[...,2],arrival)
                for time,state in ((0.,initial),(end,final)):
                    np.testing.assert_allclose(state[...,3],mass(time,arrival,v),atol=1e-14)
                errors.append(np.max(np.abs(final[...,0].ravel()-reference[:n**3,-1])))
                momentum_errors.append(np.max(np.abs(final[...,1].ravel()-reference[n**3:,-1])))
                history = np.genfromtxt(out/"history.csv",delimiter=",",names=True)
                gradient = .5*sum(np.mean((np.roll(final[...,0],-1,a)-final[...,0])**2) for a in range(3))/(box/n)**2
                expected = .5*np.mean(final[...,1]**2+final[...,3]*final[...,0]**2)+gradient
                self.assertAlmostEqual(history["total_mean"][-1],expected,places=13)
                injection = trapezoid(history["source_mean"],history["u"])
                self.assertLess(abs(history["total_mean"][-1]-history["total_mean"][0]-injection),1e-5)
                self.assertEqual(meta["scale"],1.)
            self.assertLess(errors[1],1e-5)
            self.assertGreater(errors[0]/errors[1],3.8)
            self.assertLess(momentum_errors[1],2e-5)
            self.assertGreater(momentum_errors[0]/momentum_errors[1],3.8)
            # Disabling the exact-value cache must preserve the full local state.
            uncached = root/"uncached"
            args |= {"mass_cache":"false", "result_dir":uncached}
            proc = subprocess.run([EXECUTABLE,*[f"{k}={val}" for k,val in args.items()]],capture_output=True,text=True,timeout=60)
            self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
            _, direct = read_state(uncached/"final.bin")
            np.testing.assert_array_equal(direct,final)

    def test_reject_profile_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            write_arrival_map(root/"map.afopt",np.zeros((4,4,4)),4.)
            base=dict(N=4,lSide=4,dt=.01,tMax=.1,expansion="false",arrival_map=root/"map.afopt",result_dir=root/"out")
            for change in ({"wall_width":-1},{"wall_width":"nan"},{"profile_wall_speed":0},
                           {"profile_wall_speed":1.1},{"potential":"invalid"},{"snapshot_format":"invalid"},
                           {"expansion":"true","fixedBackground":"true","energy_diagnostics":"true"}):
                proc=subprocess.run([EXECUTABLE,*[f"{k}={v}" for k,v in (base|change).items()]],capture_output=True,text=True,timeout=30)
                self.assertNotEqual(proc.returncode,0,change)
                self.assertFalse((root/"out").exists(),change)
