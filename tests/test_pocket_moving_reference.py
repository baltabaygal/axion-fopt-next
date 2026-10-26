"""Exact asymptotic wave boosts, independent of ODE mode interpolation."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from moving_packet_reference import packet_fields,left_asymptotic_fields
from planar_initial import PlanarSeed
from paper_reference import head_on_reflection


class MovingReference(unittest.TestCase):
    def test_incident_lab_phase_and_canonical_momentum(self):
        seed=PlanarSeed(speed=.3);gamma=1/np.sqrt(1-.3**2);omega=.5;k=np.sqrt(omega**2-.1**2)
        row=dict(k=k,lab_omega=omega,omega=gamma*(omega+.3*k),rest_k=gamma*(k+.3*omega),
                 reflection_real=0.,reflection_imag=0.,mass_true=1.,weight=1.,
                 x=np.array([-16.,16.]),wave=np.zeros(2,complex),derivative=np.zeros(2,complex))
        x=np.array([128.,140.,180.]);t=10.
        a,pa,grad,_=packet_fields(x,t,[row],seed)
        phase=k*(x-seed.centre)-omega*t
        np.testing.assert_allclose(a,seed.amplitude*np.cos(phase),rtol=0.,atol=2e-17)
        np.testing.assert_allclose(pa,seed.amplitude*omega*np.sin(phase),rtol=0.,atol=2e-17)
        np.testing.assert_allclose(grad,-seed.amplitude*k*np.sin(phase),rtol=0.,atol=2e-17)
        reflected=dict(row,reflection_real=1.)
        both=packet_fields(x,t,[reflected],seed)
        out=head_on_reflection(omega,.1,.3)['reflected_lab_frequency'];kout=np.sqrt(out*out-.1**2)
        phase=k*(256.-seed.centre)-kout*(x-256.)-out*t
        np.testing.assert_allclose(both[0]-a,seed.amplitude*np.cos(phase),rtol=0.,atol=5e-17)
        np.testing.assert_allclose(both[1]-pa,seed.amplitude*out*np.sin(phase),rtol=0.,atol=5e-17)

        reflected['reflected_lab_omega']=out
        asymptotic=left_asymptotic_fields(x,t,[reflected],seed)
        np.testing.assert_allclose(asymptotic[0],both[0],rtol=0.,atol=5e-17)
        np.testing.assert_allclose(asymptotic[1],both[1],rtol=0.,atol=5e-17)
