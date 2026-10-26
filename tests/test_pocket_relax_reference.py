"""Spectral force and classical mode-action checks independent of minimisation."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from model import energies
from relax_reference import RadialEnergy
from run_relaxed_pilot import field_state


class RelaxedReference(unittest.TestCase):
    def test_fixed_action_gradient_matches_energy_difference(self):
        system=RadialEnergy(.1,12.,{0:[dict(n=0,action=2.)],1:[dict(n=0,action=1.3)]})
        phi=.5*(1+np.tanh((system.r-5.)/1.2))
        direction=np.random.default_rng(8).normal(size=len(phi));direction/=np.linalg.norm(direction)
        _,gradient=system.evaluate(phi)
        estimates=[]
        for eps in (.01,.005):
            estimates.append((system.evaluate(phi+eps*direction)[0]-system.evaluate(phi-eps*direction)[0])/(2*eps))
        numerical=(4*estimates[1]-estimates[0])/3
        self.assertAlmostEqual(numerical,float(gradient@direction),delta=2e-6)

    def test_constant_mass_spectrum_matches_discrete_dirichlet_solution(self):
        h,rmax=.1,12.
        system=RadialEnergy(h,rmax,{0:[dict(n=0,action=1.)]})
        _,_,details=system.evaluate(np.ones_like(system.r),return_modes=True)
        expected=np.sqrt(4+4/h**2*np.sin(np.pi*h/(2*rmax))**2)
        self.assertAlmostEqual(details['modes'][0]['omega'],expected,places=11)

    def test_sampled_single_mode_has_assigned_action_energy(self):
        action=1.7
        system=RadialEnergy(.05,12.,{0:[dict(n=0,action=action)]})
        phi=.5*(1+np.tanh((system.r-6.)/1.))
        _,_,details=system.evaluate(phi,return_modes=True)
        modes=details['modes']
        expected=action*modes[0]['omega']
        for seed in (0,1):
            s=field_state(64,24.,system.r,phi,modes,seed)
            e=energies(s,24.,2.,.01,4.)
            actual=sum(e[k] for k in ('axion_kinetic','axion_gradient','coupling'))
            self.assertLess(abs(actual/expected-1),.01)


if __name__=='__main__':
    unittest.main()
