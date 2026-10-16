"""Optional Numba experiment checks without loading the external paper source."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np


@unittest.skipUnless(importlib.util.find_spec('numba'),'Numba required for spatial-index experiment')
class SpatialEngulfmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/stochastic'))
        import numba
        numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
        from spatial_engulfment import build_index,query_index
        cls.build=staticmethod(build_index);cls.query=staticmethod(query_index)

    def test_random_periodic_survivors_and_order(self):
        rng=np.random.default_rng(905)
        for box in (1.,6.):
            tracers=rng.uniform(-box/2,box/2,(1000,3)).astype(np.float32)
            centres=rng.uniform(-box/2,box/2,(80,3)).astype(np.float32)
            radii=rng.uniform(0,box/5,80).astype(np.float32)
            active=rng.permutation(len(tracers)).astype(np.int64)
            # Independent vectorized periodic union; random points are safely
            # away from tangency, which is checked separately below.
            delta=(tracers[:,None,:]-centres[None,:,:]).astype(np.float64)
            delta-=box*np.rint(delta/box)
            converted=np.any(np.sum(delta**2,axis=2)<=radii.astype(np.float64)[None,:]**2,axis=1)
            expected=active[~converted[active]]
            for n in (1,4,16):
                offsets,entries=self.build(centres,radii,box,n,8_000_000)
                np.testing.assert_array_equal(self.query(tracers,active,centres,radii,box,n,offsets,entries),expected)

    def test_tangency_wraparound_zero_and_large_radius(self):
        tracers=np.array([[.25,0,0],[np.nextafter(np.float32(.25),np.float32(np.inf)),0,0],
                          [-.5,-.5,-.5],[.5,.5,.5],[0,0,0]],dtype=np.float32)
        active=np.array([3,1,4,2,0],dtype=np.int64)
        for centres,radii,expected in (
            (np.array([[0,0,0]],dtype=np.float32),np.array([.25],dtype=np.float32),[3,1,2]),
            (np.array([[.5,.5,.5]],dtype=np.float32),np.array([0],dtype=np.float32),[1,4,0]),
            (np.array([[0,0,0]],dtype=np.float32),np.array([1],dtype=np.float32),[])):
            offsets,entries=self.build(centres,radii,1.,16,8_000_000)
            np.testing.assert_array_equal(self.query(tracers,active,centres,radii,1.,16,offsets,entries),np.array(expected,dtype=np.int64))

    def test_empty_and_memory_cap(self):
        tracers=np.zeros((3,3),dtype=np.float32);active=np.array([2,0],dtype=np.int64)
        centres=np.empty((0,3),dtype=np.float32);radii=np.empty(0,dtype=np.float32)
        offsets,entries=self.build(centres,radii,1.,16,8_000_000)
        np.testing.assert_array_equal(self.query(tracers,active,centres,radii,1.,16,offsets,entries),active)
        offsets,entries=self.build(tracers[:1],np.ones(1,dtype=np.float32),1.,16,0)
        self.assertEqual(len(offsets),0);self.assertEqual(len(entries),0)

    def test_first_arrival_cell_bounds_against_brute_force(self):
        from axion_fopt_next import Catalogue
        from arrival_projection import arrival_map
        rng=np.random.default_rng(981)
        for count in (0,1,100):
            cat=Catalogue(3.,.15,np.sort(rng.uniform(0,10,count)),rng.uniform(0,3,(count,3)))
            reference=cat.lattice_arrival_times(8)
            for cells in (1,4,8):
                np.testing.assert_allclose(arrival_map(cat,8,cells),reference,atol=2e-12,rtol=0)
