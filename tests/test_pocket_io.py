"""Streamed pocket binary formats and bounded-memory transfer verification."""
import sys
import struct
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
import model


class PocketIO(unittest.TestCase):
    def test_formats_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);s=np.random.default_rng(3).normal(size=(4,4,4,4))
            # Non-contiguous inputs must retain the original C-order wire format.
            s=s[::-1];source=root/'input';target=root/'state'
            model.write_initial(source,s,8.)
            payload=np.asarray(s,dtype='<f8').tobytes(order='C')
            self.assertEqual(source.read_bytes(),b'PKIN0001'+struct.pack('<Id',4,8.)+payload)
            valid=b'PKST0001'+struct.pack('<I2d',4,8.,0.)+payload
            target.write_bytes(valid)
            meta,result=model.read_state(target)
            np.testing.assert_array_equal(result,s)
            self.assertEqual(meta,dict(n=4,box=8.,u=0.))
            self.assertTrue(model.initial_transfer_exact(source,target,13))
            for bad in (valid[:-1],valid+b'x',b'BADMAGIC'+valid[8:],
                        valid[:28]+struct.pack('<d',float('nan'))+valid[36:]):
                target.write_bytes(bad)
                with self.assertRaises(ValueError):model.read_state(target)
            target.write_bytes(valid[:-1]+bytes([valid[-1]^1]))
            self.assertFalse(model.initial_transfer_exact(source,target,13))
            with self.assertRaises(ValueError):model.initial_transfer_exact(source,target,0)

    def test_prepared_and_wrapper_launch_equivalence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);s=np.zeros((4,4,4,4))
            cfg=model.config(4,8.,.01,.02,lam=2.,dv=.01,mass_tv=2.)
            launched=[]
            def fake_run(command,**kwargs):
                launched.append(Path(command[1].split('=',1)[1]).read_text())
                return type('Result',(),dict(stdout='ok',stderr='',returncode=0))()
            with patch.object(model.subprocess,'run',side_effect=fake_run):
                a=model.run(s,cfg,root/'wrapper',Path('/solver'))
                model.prepare_run(s,cfg,root/'prepared')
                b=model.run_prepared(root/'prepared',Path('/solver'))
            for key in ('initial_state','result_dir'):
                a.pop(key);b.pop(key)
            self.assertEqual(a,b)
            self.assertEqual((root/'wrapper/input.bin').read_bytes(),(root/'prepared/input.bin').read_bytes())
            self.assertEqual(len(launched),2)
            self.assertEqual((root/'prepared/run.log').read_text(),'ok')
