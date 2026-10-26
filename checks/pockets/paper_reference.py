"""Independent formation-model reference; not yet a compiled field model.

Canonical dilaton displacement phi increases from FV=0 to TV=1 in units
of its chosen field excursion. The vacuum potential is an explicit family,
not specified uniquely by the paper. Axion dependence follows PDF Eq. (2.2).
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Parameters:
    lam:float=8.
    dv:float=.01
    mass_false:float=.1
    kappa:float=4.605170185988092
    fa:float=10.

    def __post_init__(self):
        values=(self.lam,self.dv,self.mass_false,self.kappa,self.fa)
        if not np.isfinite(values).all() or self.lam<=0 or self.dv<0 or self.mass_false<=0 or self.kappa<0 or self.fa<=0:
            raise ValueError('Invalid formation-reference parameters')


def terms(phi,a,p):
    phi,a=np.asarray(phi),np.asarray(a)
    if not np.isfinite(phi).all() or not np.isfinite(a).all():
        raise ValueError('Require finite reference fields')
    with np.errstate(over='raise',invalid='raise'):
        mass2=p.mass_false**2*np.exp(p.kappa*phi)
        theta=a/p.fa
        axion=2*p.fa**2*mass2*np.sin(theta/2)**2
        vacuum=p.lam/4*phi**2*(1-phi)**2+p.dv*(1-4*phi**3+3*phi**4)
    return vacuum,axion,mass2,theta


def potential(phi,a,p):
    u,v,_,_=terms(phi,a,p);return u+v


def force(phi,a,p):
    phi=np.asarray(phi)
    _,v,m2,theta=terms(phi,a,p)
    u_phi=p.lam/2*phi*(1-phi)*(1-2*phi)+12*p.dv*phi**2*(phi-1)
    return u_phi+p.kappa*v,p.fa*m2*np.sin(theta)


def hessian(phi,a,p):
    phi=np.asarray(phi)
    _,v,m2,theta=terms(phi,a,p)
    u2=p.lam/2*(1-6*phi+6*phi**2)+12*p.dv*(3*phi**2-2*phi)
    return u2+p.kappa**2*v,p.kappa*p.fa*m2*np.sin(theta),m2*np.cos(theta)


def head_on_reflection(omega,mass_false,speed):
    """Wall moves toward incoming wave; frequencies from Appendix A kinematics.

    This is a necessary local scattering check, not an imposed reflective
    boundary or proof that all later reflected/compressed modes remain trapped.
    """
    if not np.isfinite([omega,mass_false,speed]).all() or mass_false<0 or omega<mass_false or not 0<=speed<1:
        raise ValueError('Invalid incoming wave or wall speed')
    gamma=1/np.sqrt(1-speed**2);k=np.sqrt(omega**2-mass_false**2)
    return dict(gamma=float(gamma),wall_frame_frequency=float(gamma*(omega+speed*k)),
                reflected_lab_frequency=float(gamma**2*(omega*(1+speed**2)+2*speed*k)))
