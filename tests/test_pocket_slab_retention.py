import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from slab_retention import weights

class SlabRetention(unittest.TestCase):
    def test_weight_derivative_translation_and_contraction(self):
        x=(np.arange(512)+.5)/2
        for cv,rv in ((.2,0),(0,-.1),(.2,-.1)):
            w,wd=weights(x,256.,120.,50.,cv,rv)
            eps=1e-5
            a=weights(x,256.,120.+eps*cv,50.+eps*rv,cv,rv)[0]
            b=weights(x,256.,120.-eps*cv,50.-eps*rv,cv,rv)[0]
            np.testing.assert_allclose((a-b)/(2*eps),wd,atol=1e-9,rtol=1e-7)
    def test_constant_window(self):
        w,wd=weights(np.linspace(0,256,513),256.,128.,50.,0,0)
        self.assertTrue(np.all(wd==0));self.assertTrue(np.all((w>=0)&(w<=1)))
