import sys
from pathlib import Path
import tempfile
import struct
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from field_budget import budget as original
from streaming_budget import budget


class StreamingBudget(unittest.TestCase):
    def test_matches_full_array_reference_and_slab_boundaries(self):
        rng=np.random.default_rng(94)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'state'
            for power,b in ((2,0.),(4,12.)):
                for offset in (0.,2.):
                    s=rng.normal(size=(12,12,12,4))*.4;s[...,0]+=offset
                    path.write_bytes(b'PKST0001'+struct.pack('<I2d',12,8.,.5)+s.astype('<f8').tobytes())
                    expected=original(path,lam=8.,mass_power=power,vacuum_stiffening=b,margin=.2)
                    for slab in (1,5,12):
                        observed=budget(path,lam=8.,mass_power=power,vacuum_stiffening=b,margin=.2,slab=slab)
                        self.assertEqual(observed['source_sha256'],expected['source_sha256'])
                        for key,value in expected['energies'].items():
                            self.assertAlmostEqual(observed['energies'][key],value,delta=1e-12*max(1,abs(value)))
                        for key in ('far_dilaton_energy','far_boundary_radius','instantaneous_axion_exchange'):
                            if expected[key] is None:self.assertIsNone(observed[key])
                            else:self.assertAlmostEqual(observed[key],expected[key],delta=1e-12*max(1,abs(expected[key])))

    def test_rejects_corrupt_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'state';s=np.ones((4,4,4,4));s[0,0,0,1]=np.nan
            raw=b'PKST0001'+struct.pack('<I2d',4,8.,0.)+s.astype('<f8').tobytes()
            for content in (raw,raw[:-1],raw+b'x'):
                path.write_bytes(content)
                with self.assertRaises(ValueError):budget(path)
