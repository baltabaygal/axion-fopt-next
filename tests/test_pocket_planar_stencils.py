"""Hamiltonian-force, nonlocal-current and continuum-order checks."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from planar_stencils import laplacian,gradient_density,symbol,window_budget
from planar_diagnostics import window_budget as old_budget
from paper_reference import Parameters,force


class PlanarStencils(unittest.TestCase):
    def test_fourier_symbol_positivity_and_accuracy_order(self):
        k=2*np.pi/8.;errors={2:[],4:[]}
        for n in (32,64,128):
            dx=8./n;x=(np.arange(n)+.5)*dx;f=np.sin(k*x)
            for order in (2,4):
                errors[order].append(np.max(abs(laplacian(f,dx,order)+k*k*f)))
                np.testing.assert_allclose(laplacian(f,dx,order),-symbol(k,dx,order)*f,rtol=0.,atol=2e-12)
                self.assertTrue(np.all(symbol(np.linspace(0,np.pi/dx,256),dx,order)>=0.))
        for order,expected in ((2,4.),(4,16.)):
            self.assertGreater(errors[order][0]/errors[order][1],expected*.95)
        self.assertAlmostEqual(symbol(np.pi/.25,.25,4),16/(3*.25**2))

    def test_gradient_hamiltonian_derivative_is_negative_laplacian(self):
        rng=np.random.default_rng(123);f=rng.normal(size=32);dx=.3;eps=1e-6
        for order in (2,4):
            self.assertGreater(np.sum(gradient_density(f,dx,order)),0.)
            for j in (0,4,31):
                plus=f.copy();minus=f.copy();plus[j]+=eps;minus[j]-=eps
                derivative=(np.sum(gradient_density(plus,dx,order))-np.sum(gradient_density(minus,dx,order)))/(2*eps)
                self.assertAlmostEqual(derivative,-laplacian(f,dx,order)[j],delta=2e-7)

    def test_signed_long_bond_currents_close_both_sector_balances(self):
        n=32;box=8.;dx=box/n;x=(np.arange(n)+.5)*dx;k=2*np.pi/box
        s=np.column_stack((.5+.3*np.sin(k*x),.3*np.cos(k*x),.05*np.cos(k*x)+.02,-.1*np.sin(k*x)))
        p=Parameters();fp,fa=force(s[:,0],s[:,1],p)
        w=.5+.3*np.sin(k*x);wd=-.2*.3*k*np.cos(k*x);eps=1e-6
        for order in (2,4):
            rhs=np.column_stack((s[:,2],s[:,3],laplacian(s[:,0],dx,order)-fp,laplacian(s[:,1],dx,order)-fa))
            v=window_budget(s,box,w,wd,p,order)
            plus=window_budget(s+eps*rhs,box,w+eps*wd,wd,p,order)
            minus=window_budget(s-eps*rhs,box,w-eps*wd,wd,p,order)
            for key,expected in (('energy',-v['phi_flux']-v['axion_flux']),('phi_energy',-v['phi_flux']-v['exchange']),('axion_energy',-v['axion_flux']+v['exchange'])):
                self.assertAlmostEqual((plus[key]-minus[key])/(2*eps),expected,delta=3e-9)
            if order==2:
                old=old_budget(s,box,w,wd,p)
                for key in old:self.assertAlmostEqual(v[key],old[key],delta=1e-14)
