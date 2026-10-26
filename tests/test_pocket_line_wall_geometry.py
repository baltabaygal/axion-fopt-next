"""Actual C++ wall diagnostic vs independent Python roots, including the cut."""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks/pockets'))
from planar_initial import PlanarSeed,prepare
from planar_diagnostics import wall_geometry


@unittest.skipUnless(shutil.which('c++'),'C++ compiler unavailable')
class LineWallGeometry(unittest.TestCase):
    def test_cpp_roots_velocity_width_and_ambiguity(self):
        code=r'''
#include "line_wall_geometry.h"
#include <iostream>
#include <iomanip>
int main(){int n;double box;std::cin>>n>>box;std::vector<double>s(4*n);
for(double&x:s)std::cin>>x;
try{auto v=lineWallGeometry(s,n,box);std::cout<<std::setprecision(17);for(double x:v)std::cout<<x<<" ";}
catch(const std::runtime_error&){std::cout<<"ambiguous";}}
'''
        with tempfile.TemporaryDirectory(prefix='line-geometry-') as folder:
            path=Path(folder);(path/'probe.cpp').write_text(code)
            build=subprocess.run(['c++','-std=c++17','-O0','-I',str(ROOT/'checks/pockets'),str(path/'probe.cpp'),'-o',str(path/'probe')],capture_output=True,text=True,timeout=20)
            self.assertEqual(build.returncode,0,build.stderr)
            def probe(s,box):
                text=f'{len(s)} {box}\n'+' '.join(format(v,'.17g') for v in s.ravel())
                return subprocess.run([str(path/'probe')],input=text,capture_output=True,text=True,check=True,timeout=5).stdout
            for speed in (0.,.3):
                _,s=prepare(2048,PlanarSeed(speed=speed))
                for shift in (0,1024):
                    shifted=np.roll(s,shift,axis=0);expected=wall_geometry(shifted,512.)
                    actual=np.fromstring(probe(shifted,512.),sep=' ')
                    np.testing.assert_allclose(actual,list(expected.values()),rtol=0.,atol=1e-13)
            bad=np.zeros((64,4));bad[:,0]=.5+.4*np.sin(np.arange(64)*6*np.pi/64)
            self.assertEqual(probe(bad,8.),'ambiguous')
