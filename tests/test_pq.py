"""PQ force/topology contracts and optional compiled-backend reference validation."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
try:
    import scipy
except ImportError:
    scipy=None
if scipy is not None:
    spec=importlib.util.spec_from_file_location('pq_reference',ROOT/'checks/pq_strings/model.py')
    pq=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pq)

@unittest.skipIf(scipy is None,'PQ validation requires optional scipy')
class PQContracts(unittest.TestCase):
    def test_force_origin_and_vacuum(self):
        self.assertTrue(np.array_equal(pq.force(np.array(0.),np.array(0.),.4,.03),[0.,0.]))
        direction=np.linspace(-np.pi,np.pi,31)
        gradient=pq.force(1e-12*np.cos(direction),1e-12*np.sin(direction),.4,.03)
        self.assertLess(np.max(np.linalg.norm(gradient,axis=-1)),1e-12)
        self.assertTrue(np.array_equal(pq.force(np.array(1.),np.array(0.),.4,.03),[0.,0.]))

    def test_force_is_potential_gradient(self):
        rng=np.random.default_rng(19)
        x,y=rng.uniform(-2,2,(2,100));epsilon=1e-6
        difference=np.stack(((pq.potential(x+epsilon,y,.4,.03)-pq.potential(x-epsilon,y,.4,.03))/(2*epsilon),
                             (pq.potential(x,y+epsilon,.4,.03)-pq.potential(x,y-epsilon,.4,.03))/(2*epsilon)),axis=-1)
        np.testing.assert_allclose(pq.force(x,y,.4,.03),difference,atol=2e-8,rtol=2e-8)

    def test_periodic_neutral_strings_and_no_wall_at_uniform_vacuum(self):
        state=pq.deterministic(16,8.)
        d,w=pq.defects(state,8.,20.)
        self.assertEqual(d['plaquettes'],4*16)
        self.assertEqual(d['signed_winding'],[0,0,0])
        self.assertEqual(np.count_nonzero(w[0][...,0]),4)
        uniform=np.zeros_like(state);uniform[...,0]=1.
        d,_=pq.defects(uniform,8.,20.)
        self.assertEqual(d['plaquettes'],0)
        self.assertEqual(d['wall_links'],0)
        phase=2*np.pi*(np.arange(16)+.5)/16
        sheet=np.zeros_like(state)
        sheet[...,0]=np.cos(phase)[None,None,:]
        sheet[...,1]=np.sin(phase)[None,None,:]
        d,_=pq.defects(sheet,8.,20.)
        self.assertEqual(d['plaquettes'],0)
        self.assertEqual(d['wall_links'],16**2)
        with self.assertRaises(ValueError):
            pq.defects(np.zeros_like(state),8.,20.)

    def test_mass_hubble_equality(self):
        cfg=pq.config(4,4.,.01,1.,origin=20.,star=27.,exponent=8.)
        tau,a,chi=pq.coefficients(cfg,7.)
        self.assertAlmostEqual(np.sqrt(chi),20/27**2)
        self.assertAlmostEqual(a,27/20)

    @unittest.skipUnless(os.environ.get('AFOPT_PQ_EXECUTABLE'),'Set AFOPT_PQ_EXECUTABLE')
    def test_compiled_evolved_wave_yield(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([sys.executable,str(ROOT/'checks/pq_strings/validate_yield.py'),
                '--output',str(Path(tmp)/'validation'),'--executable',os.environ['AFOPT_PQ_EXECUTABLE']],
                capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    @unittest.skipUnless(os.environ.get('AFOPT_PQ_EXECUTABLE'),'Set AFOPT_PQ_EXECUTABLE')
    def test_compiled_core_growth_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([sys.executable,str(ROOT/'checks/pq_strings/validate_preparation.py'),
                '--output',str(Path(tmp)/'validation'),'--executable',os.environ['AFOPT_PQ_EXECUTABLE']],
                capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    @unittest.skipUnless(os.environ.get('AFOPT_PQ_EXECUTABLE'),'Set AFOPT_PQ_EXECUTABLE')
    def test_compiled_reference_suite(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'validation'
            result=subprocess.run([sys.executable,str(ROOT/'checks/pq_strings/validate.py'),
                                   '--output',str(output),'--executable',os.environ['AFOPT_PQ_EXECUTABLE']],
                                  capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            report=json.loads((output/'validation.json').read_text())
            self.assertTrue(report['passed'])

    @unittest.skipUnless(os.environ.get('AFOPT_PQ_EXECUTABLE'),'Set AFOPT_PQ_EXECUTABLE')
    def test_compiled_capped_mass_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([sys.executable,str(ROOT/'checks/pq_strings/validate_cap.py'),
                '--output',str(Path(tmp)/'validation'),'--executable',os.environ['AFOPT_PQ_EXECUTABLE']],
                capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
