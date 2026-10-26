import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from wall_kinematics import measure,summary
from wall_diagnostics import local_widths


class WallKinematics(unittest.TestCase):
    def test_translating_contracted_planar_profiles(self):
        dx=.025;x=(np.arange(512)+.5)*dx-6.4;w=.5
        for speed in (0.,.6,.85):
            gamma=1/np.sqrt(1-speed**2)
            f=.5*(1+np.tanh(gamma*(abs(x)-3)/w))
            derivative=2*gamma/w*f*(1-f)*np.sign(x)
            phi=np.repeat(f[:,None],8,axis=1);pi=np.repeat((-speed*derivative)[:,None],8,axis=1)
            wall,data=measure(phi,pi,dx,reach=2.)
            self.assertEqual(wall,local_widths(phi,dx,reach=2.))
            result=summary(data,dx)
            self.assertAlmostEqual(result['all']['median_abs_level_speed'],speed,delta=.01)
            self.assertAlmostEqual(result['median_boost_corrected_width'],2*np.arctanh(.8)*w,delta=.03)

    def test_superluminal_level_speeds_are_retained(self):
        r=summary(np.array([[2.,1.3,0.,0.],[6.,.2,1.,1.]]),.1)
        self.assertEqual(r['thin']['fraction_abs_speed_ge_one'],1.)
        self.assertEqual(r['minimum_width_abs_level_speed'],1.3)
