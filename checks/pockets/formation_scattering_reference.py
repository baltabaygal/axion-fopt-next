"""Independent continuum ODE scattering of a linear axion across a planar wall.

Weak-field limit only: phi is the static degenerate vacuum kink, not an
axion-driven wall. This is not a lattice evolution or a survival experiment.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp


def scatter(omega,mass_false=.1,kappa=math.log(100.),width=1.,extent=16.,
            rtol=2e-11,atol=2e-13,profile='tanh',samples=257):
    values=[omega,mass_false,kappa,width,extent,rtol,atol]
    if (not np.isfinite(values).all() or mass_false<=0 or kappa<0 or width<=0 or
        extent<=0 or omega<=mass_false or rtol<=0 or atol<=0 or
        profile not in ('tanh','step') or samples<2):
        raise ValueError('Invalid planar scattering parameters')
    mass_true=mass_false*math.exp(kappa/2.)
    if not math.isfinite(mass_true):raise ValueError('Nonfinite exterior mass')
    left_k=math.sqrt(omega**2-mass_false**2)
    # Outgoing wave above threshold; decaying tail below it. exp(i k x)
    # with positive imaginary k gives the correct evanescent sign.
    right_k=complex(np.sqrt(complex(omega**2-mass_true**2)))
    def rhs(x,y):
        phi=.5*(1+math.tanh(x/width)) if profile=='tanh' else float(x>=0.)
        mass2=mass_false**2*math.exp(kappa*phi)
        return [y[1],(mass2-omega**2)*y[0]]
    x=np.linspace(extent,-extent,samples)
    sol=solve_ivp(rhs,(extent,-extent),[1.+0j,1j*right_k],method='DOP853',
                  rtol=rtol,atol=atol,t_eval=x)
    if not sol.success:raise RuntimeError(sol.message)
    y,dy=sol.y[:,-1];left=-extent
    incoming=.5*(y+dy/(1j*left_k))*np.exp(-1j*left_k*left)
    outgoing=.5*(y-dy/(1j*left_k))*np.exp(1j*left_k*left)
    reflection=outgoing/incoming
    reflected_flux=float(abs(reflection)**2)
    transmitted_flux=float(right_k.real/(left_k*abs(incoming)**2))
    return dict(omega=omega,mass_false=mass_false,mass_true=mass_true,width=width,
        extent=extent,rtol=rtol,atol=atol,profile=profile,nfev=sol.nfev,
        reflection_real=float(reflection.real),reflection_imag=float(reflection.imag),
        reflected_flux=reflected_flux,transmitted_flux=transmitted_flux,
        flux_residual=abs(reflected_flux+transmitted_flux-1.),
        x=sol.t[::-1],wave=(sol.y[0]/incoming)[::-1],derivative=(sol.y[1]/incoming)[::-1])


def step_amplitudes(omega,mass_false,mass_true):
    """Continuity of psi and psi' at a zero-width real mass step."""
    kl=np.sqrt(complex(omega**2-mass_false**2));kr=np.sqrt(complex(omega**2-mass_true**2))
    if kl.real<=0:raise ValueError('Incoming propagating mode required')
    r=(kl-kr)/(kl+kr);t=2*kl/(kl+kr)
    return r,float(abs(r)**2),float(kr.real/kl.real*abs(t)**2)
