"""Formation-model semi-discrete fixed and moving sector continuity."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from model import radius,laplacian
from paper_reference import Parameters,force
from monitor_reference import monitor
from tracking_monitor_reference import tracking_monitor


class FormationMonitors(unittest.TestCase):
    def test_fixed_and_moving_sector_exchange_and_current(self):
        n=12;box=8.;dx=box/n;r=radius(n,box)
        state=np.zeros((n,n,n,4))
        state[...,0]=.5*(1+np.tanh((r-1.7)/.6))
        state[...,1]=2.*np.exp(-r*r/3.)
        state[...,2]=-.1*(np.roll(state[...,0],-1,0)-np.roll(state[...,0],1,0))
        state[...,2]+=.02*np.exp(-r*r/3.)
        state[...,3]=.15*np.exp(-r*r/3.)
        p=Parameters();fp,fa=force(state[...,0],state[...,1],p)
        rhs=np.stack((state[...,2],state[...,3],laplacian(state[...,0],dx)-fp,
                      laplacian(state[...,1],dx)-fa),axis=-1)
        for prefix,measure in (
            ('monitor',lambda s:monitor(s,box,2.7,p.lam,p.dv,1.,formation=p)),
            ('tracking',lambda s:tracking_monitor(s,box,2.4,3.5,p.lam,p.dv,1.,formation=p))):
            eps=1e-6;v=measure(state);plus=measure(state+eps*rhs);minus=measure(state-eps*rhs)
            self.assertGreater(abs(v[prefix+'_exchange']),1e-5)
            flux_phi=v[prefix+'_phi_flux'];flux_a=v[prefix+'_axion_flux'];q=v[prefix+'_exchange']
            for sector,expected in (('energy',-flux_phi-flux_a),
                                    ('phi_energy',-flux_phi-q),('axion_energy',-flux_a+q)):
                observed=(plus[prefix+'_'+sector]-minus[prefix+'_'+sector])/(2*eps)
                self.assertLess(abs(observed-expected)/v[prefix+'_energy'],1e-7)
            if prefix=='tracking':self.assertGreater(np.linalg.norm(v['velocity']),.01)
