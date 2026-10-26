import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from compare_field_grids import compare


class GridComparison(unittest.TestCase):
    @staticmethod
    def history(times):
        names=('u','fv_radius','axion_far','axion_kinetic','axion_gradient','coupling','monitor_energy',
               'tracking_energy','tracking_phi_flux','tracking_axion_flux','tracking_exchange',
               'tracking_cx','tracking_cy','tracking_cz')
        h=np.zeros(len(times),dtype=[(name,float) for name in names]);h['u']=times
        h['fv_radius']=20.;h['axion_kinetic']=1-.01*h['u'];h['axion_gradient']=1.
        h['monitor_energy']=2-.003*h['u'];h['tracking_energy']=2-.003*h['u']
        h['tracking_phi_flux']=.002;h['tracking_axion_flux']=.001;h['tracking_exchange']=.0005
        h['tracking_cx']=63.9
        return h

    def test_prefix_uses_shared_endpoint_and_periodic_centroids(self):
        a=self.history([0,1,2]);b=self.history([0,.5,1]);b['tracking_cx']=.1
        r=compare(a,b,64.)
        self.assertEqual(r['depletion_difference'],0.)
        self.assertEqual(r['monitor_retention']['coarse'],r['monitor_retention']['fine'])
        self.assertFalse(r['history_gates']['full_duration'])
        self.assertTrue(all(r['moving_comparison']['gates'].values()))
        self.assertAlmostEqual(r['moving_comparison']['max_centroid_distance'],.2)

    def test_detects_storage_and_current_discrepancies(self):
        a=self.history([0,40,80]);b=a.copy()
        b['tracking_energy'][-1]+=.01;b['tracking_phi_flux']+=.0001
        b['tracking_cx']=.3
        r=compare(a,b,64.)
        self.assertTrue(r['history_gates']['full_duration'])
        for gate in ('storage','current','centroid'):self.assertFalse(r['moving_comparison']['gates'][gate])

    def test_rejects_duplicate_or_missing_initial_times(self):
        a=self.history([0,1])
        for times in ([0,0],[1,2]):
            with self.assertRaises(ValueError):compare(a,self.history(times),64.)
