"""Finite-difference energy continuity checks for the independent monitor."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from model import force,laplacian
from monitor_reference import monitor,continuity_summary


class MonitorContinuity(unittest.TestCase):
    def test_moving_summary_uses_corrected_current_and_sector_exchange(self):
        names=('energy','phi_energy','axion_energy','phi_flux','axion_flux','exchange')
        h=np.zeros(9,dtype=[('u',float)]+[('tracking_'+name,float) for name in names])
        t=np.linspace(0,2,9);h['u']=t
        h['tracking_phi_energy']=60-5*t;h['tracking_axion_energy']=40+t
        h['tracking_energy']=100-4*t
        h['tracking_phi_flux']=3;h['tracking_axion_flux']=1;h['tracking_exchange']=2
        result=continuity_summary(h,'tracking')
        self.assertEqual(result['max_continuity_residual'],dict(total=0.,phi=0.,axion=0.))
        self.assertEqual(result['net_outward_phi_energy'],6.)
        self.assertEqual(result['net_energy_into_axions'],4.)
        self.assertIsNone(continuity_summary(h))

    def test_sector_continuity_including_coupling_exchange(self):
        rng=np.random.default_rng(19);n=8;box=8.
        state=rng.normal(size=(n,n,n,4))*.15
        state[...,0]+=.5
        for power in (2,4):
            fphi,fa=force(state[...,0],state[...,1],2.,.01,4.,power)
            rhs=np.stack((state[...,2],state[...,3],laplacian(state[...,0],1.)-fphi,
                          laplacian(state[...,1],1.)-fa),-1)
            eps=1e-6
            plus=monitor(state+eps*rhs,box,2.7,2.,.01,4.,power)
            minus=monitor(state-eps*rhs,box,2.7,2.,.01,4.,power)
            reference=monitor(state,box,2.7,2.,.01,4.,power)
            derivative={k:(plus[k]-minus[k])/(2*eps) for k in ('monitor_energy','monitor_phi_energy','monitor_axion_energy')}
            self.assertAlmostEqual(derivative['monitor_phi_energy'],
                -reference['monitor_phi_flux']-reference['monitor_exchange'],delta=1e-7)
            self.assertAlmostEqual(derivative['monitor_axion_energy'],
                -reference['monitor_axion_flux']+reference['monitor_exchange'],delta=1e-7)
            self.assertAlmostEqual(derivative['monitor_energy'],
                -reference['monitor_phi_flux']-reference['monitor_axion_flux'],delta=1e-7)


if __name__=='__main__':
    unittest.main()
