import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checks/pockets'))
from paper_reference import Parameters,potential,force,hessian,head_on_reflection


class PaperReference(unittest.TestCase):
    def test_force_and_hessian_by_independent_differences(self):
        p=Parameters();eps=1e-5
        for x,a in ((0.,.4),(.3,1.2),(.8,-2.),(1.,.7)):
            f,g=force(x,a,p);xx,xy,aa=hessian(x,a,p)
            numerical=((potential(x+eps,a,p)-potential(x-eps,a,p))/(2*eps),
                       (potential(x,a+eps,p)-potential(x,a-eps,p))/(2*eps))
            np.testing.assert_allclose([f,g],numerical,rtol=1e-7,atol=1e-8)
            np.testing.assert_allclose([xx,xy,aa],[(force(x+eps,a,p)[0]-force(x-eps,a,p)[0])/(2*eps),
                (force(x,a+eps,p)[0]-force(x,a-eps,p)[0])/(2*eps),
                (force(x,a+eps,p)[1]-force(x,a-eps,p)[1])/(2*eps)],rtol=1e-7,atol=1e-8)

    def test_vacua_periodicity_and_false_vacuum_backreaction(self):
        p=Parameters();self.assertAlmostEqual(potential(0.,0.,p),p.dv)
        self.assertAlmostEqual(potential(1.,0.,p),0.)
        self.assertAlmostEqual(hessian(1.,0.,p)[2]/hessian(0.,0.,p)[2],100.)
        self.assertAlmostEqual(potential(.4,.3,p),potential(.4,.3+2*np.pi*p.fa,p))
        self.assertGreater(force(0.,.5,p)[0],0.)

    def test_quadratic_limit_and_no_exponential_clipping(self):
        p=Parameters(fa=1e5);a=.3
        expected=.5*p.mass_false**2*np.exp(p.kappa*.4)*a*a
        self.assertAlmostEqual(potential(.4,a,p)-potential(.4,0.,p),expected,delta=1e-12)
        with self.assertRaises(FloatingPointError):potential(1000.,a,p)

    def test_head_on_massless_doppler_limits(self):
        r=head_on_reflection(.6,0.,.6)
        self.assertAlmostEqual(r['wall_frame_frequency'],1.2)
        self.assertAlmostEqual(r['reflected_lab_frequency'],2.4)
        r=head_on_reflection(.6,0.,.999)
        self.assertGreater(r['wall_frame_frequency'],26.)
        self.assertAlmostEqual(r['reflected_lab_frequency'],.6*1999,places=7)

    def test_array_fields_and_nonfinite_rejection(self):
        p=Parameters();f,g=force([0.,.3],[.1,.2],p)
        self.assertEqual(f.shape,(2,));self.assertEqual(g.shape,(2,))
        with self.assertRaises(ValueError):potential(np.nan,0.,p)
        with self.assertRaises(ValueError):force(0.,np.inf,p)
