"""Small scalar-kernel cross-check; no framework build or lattice evolution."""
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks/pockets'))
from paper_reference import Parameters,terms,force,hessian


@unittest.skipUnless(shutil.which('c++'),'C++ compiler unavailable')
class FormationCpp(unittest.TestCase):
    def test_compiled_kernel_matches_independent_reference(self):
        source=r'''
#include "pocket_potential.h"
#include <iostream>
#include <iomanip>
int main() {
  pockets::FormationParameters p;
  double phi,a;
  std::cout << std::setprecision(17);
  while(std::cin>>phi>>a) {
    try {
      const auto v=pockets::formationPotential(phi,a,p);
      if(TempLat::pocketValue<0>(phi,a,8.,.01,1.,2,0.,&p)!=v.potential() ||
         TempLat::pocketValue<1>(phi,a,8.,.01,1.,2,0.,&p)!=v.dphi ||
         TempLat::pocketValue<2>(phi,a,8.,.01,1.,2,0.,&p)!=v.da)return 4;
      std::cout<<v.vacuum<<" "<<v.axion<<" "<<v.mass_squared<<" "
        <<v.dphi<<" "<<v.da<<" "<<v.hphi_phi<<" "<<v.hphi_a<<" "
        <<v.ha_a<<" "<<v.axion_exchange(.37)<<"\n";
    } catch(const std::overflow_error&) { std::cout<<"overflow\n"; }
  }
  p.fa=0.;
  try { pockets::formationPotential(0.,0.,p); return 2; }
  catch(const std::invalid_argument&) {}
}
'''
        rng=np.random.default_rng(240701676)
        points=[(0.,0.),(1.,0.),(0.,.5),(1.-1e-7,.3),(.4,2*np.pi*10+.3)]
        points+=list(zip(rng.uniform(-.2,1.2,64),rng.uniform(-40,40,64)))
        with tempfile.TemporaryDirectory(prefix='formation-kernel-') as directory:
            directory=Path(directory)
            (directory/'probe.cpp').write_text(source)
            compilation=subprocess.run(['c++','-std=c++17','-O0','-I',str(ROOT/'checks/pockets'),
                            str(directory/'probe.cpp'),'-o',str(directory/'probe')],
                           capture_output=True,text=True,timeout=30)
            self.assertEqual(compilation.returncode,0,compilation.stderr)
            data=''.join(f'{x:.17g} {a:.17g}\n' for x,a in points)+'1000 .3\n'
            result=subprocess.run([str(directory/'probe')],input=data,capture_output=True,
                                  text=True,check=True,timeout=10).stdout.splitlines()
        self.assertEqual(result[-1],'overflow')
        p=Parameters()
        for (x,a),line in zip(points,result[:-1]):
            u,v,m2,_=terms(x,a,p)
            # Reference U is expanded and loses relative precision near phi=1.
            expected=[u,v,m2,*force(x,a,p),*hessian(x,a,p),p.kappa*v*.37]
            np.testing.assert_allclose(np.fromstring(line,sep=' '),expected,rtol=2e-12,atol=2e-14)
        self.assertEqual(len(result),len(points)+1)
