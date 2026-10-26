import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from model import radius,energies
from tracking_monitor_reference import tracking_monitor
from window_support_audit import audit


class WindowSupport(unittest.TestCase):
    def test_matches_weighted_monitor_and_global_energy(self):
        r=radius(16,8.,[.6,4.2,4.1]);state=np.zeros((16,16,16,4))
        state[...,0]=.5*(1+np.tanh((r-1.5)/.4))
        state[...,1]=.2*np.exp(-r*r/2);state[...,3]=.1*np.exp(-r*r/3)
        for b in (0.,12.):
            ref=tracking_monitor(state,8.,2.5,3.5,2.,.01,4.,4,b)
            for slab in (1,5,16):
                result=audit(state,8.,ref['centre'],2.5,3.5,2.,.01,4.,4,b,slab)
                self.assertAlmostEqual(result['weighted_total'],ref['tracking_energy'],delta=1e-12)
                total=energies(state,8.,2.,.01,4.,mass_power=4,vacuum_stiffening=b)['total']
                self.assertAlmostEqual(sum(result[k+'_total'] for k in ('inner','annulus','outside')),total,delta=1e-12)
                self.assertLess(result['threshold_max_radii']['0.5'],2.5)
