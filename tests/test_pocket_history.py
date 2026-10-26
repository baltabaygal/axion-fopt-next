"""Check tensor assembly using independent ellipsoids and rotations."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from history_diagnostics import asphericity


class PocketHistory(unittest.TestCase):
    def test_ellipsoid_and_rotated_ellipsoid_have_same_shape(self):
        dtype=[(key,float) for key in ('sxx','syy','szz','sxy','sxz','syz','fv_volume')]
        h=np.zeros(3,dtype=dtype);h['fv_volume']=1.
        angle=.71
        rotation=np.array([[np.cos(angle),-np.sin(angle),0.],
                           [np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])
        covariance=np.diag([1.,2.,6.])
        for i,matrix in enumerate((covariance,rotation@covariance@rotation.T,np.eye(3)*3)):
            for key,(a,b) in zip(('sxx','syy','szz','sxy','sxz','syz'),
                                 ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))):
                h[key][i]=matrix[a,b]
        np.testing.assert_allclose(asphericity(h),[5/3,5/3,0.],atol=1e-14)
        h['fv_volume'][1]=0
        self.assertTrue(np.isnan(asphericity(h)[1]))


if __name__=='__main__':
    unittest.main()
