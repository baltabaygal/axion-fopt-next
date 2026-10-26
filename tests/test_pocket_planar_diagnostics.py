import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from planar_diagnostics import wall_geometry,window_budget
from planar_initial import prepare,PlanarSeed
from paper_reference import Parameters,force


class PlanarDiagnostics(unittest.TestCase):
    def test_boosted_wall_geometry_and_ambiguity(self):
        for speed in (0.,.3):
            _,state=prepare(4096,PlanarSeed(speed=speed))
            result=wall_geometry(state,512.)
            self.assertAlmostEqual(result['position'],256.,delta=1e-12)
            self.assertAlmostEqual(result['phase_velocity'],-speed,delta=1e-3)
            self.assertAlmostEqual(result['width_10_90'],2*np.arctanh(.8)*np.sqrt(1-speed**2),delta=.01)
        bad=np.zeros((64,4));bad[:,0]=.5+.4*np.sin(np.arange(64)*6*np.pi/64)
        with self.assertRaises(ValueError):wall_geometry(bad,8.)

    def test_exact_moving_window_sector_identity(self):
        n=32;box=8.;dx=box/n;x=(np.arange(n)+.5)*dx;k=2*np.pi/box
        s=np.column_stack((.5+.4*np.sin(k*x),.3*np.cos(k*x),.05*np.cos(k*x)+.02,-.1*np.sin(k*x)))
        p=Parameters();fp,fa=force(s[:,0],s[:,1],p)
        lap=(np.roll(s[:,:2],-1,axis=0)+np.roll(s[:,:2],1,axis=0)-2*s[:,:2])/dx**2
        rhs=np.column_stack((s[:,2:],lap-np.column_stack((fp,fa))))
        w=.5+.3*np.sin(k*x);wd=-.2*.3*k*np.cos(k*x);eps=1e-6
        value=window_budget(s,box,w,wd,p)
        plus=window_budget(s+eps*rhs,box,w+eps*wd,wd,p)
        minus=window_budget(s-eps*rhs,box,w-eps*wd,wd,p)
        self.assertGreater(abs(value['exchange']),1e-5)
        for key,expected in (('phi_energy',-value['phi_flux']-value['exchange']),
                             ('axion_energy',-value['axion_flux']+value['exchange']),
                             ('energy',-value['phi_flux']-value['axion_flux'])):
            self.assertAlmostEqual((plus[key]-minus[key])/(2*eps),expected,delta=2e-9)
