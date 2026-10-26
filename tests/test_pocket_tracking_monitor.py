"""Moving-volume continuity and periodic translation identities."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from model import force,laplacian,radius
from tracking_monitor_reference import tracking_monitor


class TrackingMonitor(unittest.TestCase):
    @staticmethod
    def state():
        n=16;r=radius(n,8.);s=np.zeros((n,n,n,4))
        s[...,0]=.5*(1+np.tanh((r-1.8)/.6))
        s[...,1]=.15*np.exp(-r*r/2)
        s[...,3]=.08*np.exp(-r*r/3)
        s[...,2]=-.12*(np.roll(s[...,0],-1,0)-np.roll(s[...,0],1,0))
        s[...,2]+=.02*np.random.default_rng(22).normal(size=(n,n,n))
        return s

    def test_exact_sector_identity_in_a_moving_window(self):
        s=self.state();eps=1e-6
        for power in (2,4):
            for stiffening in (0.,12.):
                fphi,fa=force(s[...,0],s[...,1],2.,.01,4.,power,stiffening)
                rhs=np.stack((s[...,2],s[...,3],laplacian(s[...,0],.5)-fphi,
                              laplacian(s[...,1],.5)-fa),-1)
                parameters=(8.,2.5,3.5,2.,.01,4.,power,stiffening)
                initial=tracking_monitor(s,*parameters)
                plus=tracking_monitor(s+eps*rhs,*parameters);minus=tracking_monitor(s-eps*rhs,*parameters)
                expected=dict(tracking_phi_energy=-initial['tracking_phi_flux']-initial['tracking_exchange'],
                              tracking_axion_energy=-initial['tracking_axion_flux']+initial['tracking_exchange'],
                              tracking_energy=-initial['tracking_phi_flux']-initial['tracking_axion_flux'])
                self.assertGreater(np.linalg.norm(initial['velocity']),.01)
                for key,value in expected.items():
                    observed=(plus[key]-minus[key])/(2*eps)
                    self.assertLess(abs(observed-value)/initial['tracking_energy'],1e-7)

    def test_periodic_translation_preserves_energy_current_and_velocity(self):
        state=self.state();shifted=np.roll(state,(7,-6,5),axis=(0,1,2))
        a=tracking_monitor(state,8.,2.5,3.5,2.,.01,4.)
        b=tracking_monitor(shifted,8.,2.5,3.5,2.,.01,4.)
        for key in ('tracking_energy','tracking_phi_energy','tracking_axion_energy',
                    'tracking_phi_flux','tracking_axion_flux','tracking_exchange'):
            self.assertAlmostEqual(a[key],b[key],delta=1e-12)
        np.testing.assert_allclose(a['velocity'],b['velocity'],atol=1e-12)


if __name__=='__main__':
    unittest.main()
