"""Moving-window semi-discrete current reference; no field evolution."""
import numpy as np

from model import dilaton_potential
from pocket_geometry import periodic_centroid


def tracking_monitor(state,box,inner,outer,lam,dv,mass2,power=4,vacuum_stiffening=0.,formation=None):
    if not 0<inner<outer<box/2:
        raise ValueError('Require 0<inner<outer<L/2')
    n=len(state);dx=box/n;phi,a,pphi,pa=(state[...,j] for j in range(4))
    phase=np.clip((phi-.1)/.8,0.,1.)
    eta=1-phase**3*(10+phase*(-15+6*phase))
    eta_dot=-30*phase**2*(1-phase)**2/.8*pphi
    centre,coherence=periodic_centroid(eta,box);weight=eta.sum()
    x=(np.arange(n)+.5)*dx;coordinates=[];velocity=[]
    for axis in range(3):
        d=(x-centre[axis]+box/2)%box-box/2;coordinates.append(d)
        moments=eta_dot.sum(axis=tuple(j for j in range(3) if j!=axis))
        velocity.append(float(moments@d/weight))
    X,Y,Z=coordinates[0][:,None,None],coordinates[1][None,:,None],coordinates[2][None,None,:]
    r=np.sqrt(X*X+Y*Y+Z*Z);s=np.clip((r-inner)/(outer-inner),0.,1.)
    window=1-s**3*(10+s*(-15+6*s))
    radial_derivative=-30*s*s*(1-s)**2/(outer-inner)
    numerator=velocity[0]*X+velocity[1]*Y+velocity[2]*Z
    window_dot=-radial_derivative*np.divide(numerator,r,out=np.zeros_like(r),where=r>0)
    phi_density=.5*pphi**2+dilaton_potential(phi,lam,dv,vacuum_stiffening)
    a_density=.5*pa**2+.5*mass2*phi**power*a*a
    if formation is not None:
        from paper_reference import terms
        if formation.lam!=lam or formation.dv!=dv or vacuum_stiffening!=0.:
            raise ValueError('Inconsistent formation vacuum parameters')
        u,va,_,_=terms(phi,a,formation)
        phi_density=.5*pphi**2+u;a_density=.5*pa**2+va
    fluxes=[]
    for field,momentum,density in ((phi,pphi,phi_density),(a,pa,a_density)):
        flux=0.
        for axis in range(3):
            following=np.roll(field,-1,axis);previous=np.roll(field,1,axis)
            density+=((following-field)**2+(field-previous)**2)/(4*dx*dx)
            current=-(following-field)*(momentum+np.roll(momentum,-1,axis))/(2*dx)
            flux+=float(np.sum(current*(window-np.roll(window,-1,axis)))*dx*dx)
        fluxes.append(flux)
    ephi=float(np.sum(window*phi_density)*dx**3);ea=float(np.sum(window*a_density)*dx**3)
    adv_phi=float(np.sum(window_dot*phi_density)*dx**3);adv_a=float(np.sum(window_dot*a_density)*dx**3)
    exchange=float(np.sum(window*.5*power*mass2*phi**(power-1)*a*a*pphi)*dx**3)
    if formation is not None:
        exchange=float(np.sum(window*formation.kappa*va*pphi)*dx**3)
    return dict(tracking_energy=ephi+ea,tracking_phi_energy=ephi,tracking_axion_energy=ea,
                tracking_phi_flux=fluxes[0]-adv_phi,tracking_axion_flux=fluxes[1]-adv_a,
                tracking_exchange=exchange,tracking_lab_phi_flux=fluxes[0],tracking_lab_axion_flux=fluxes[1],
                tracking_phi_advection=adv_phi,tracking_axion_advection=adv_a,
                centre=centre.tolist(),velocity=velocity,centroid_coherence=coherence,
                phase_volume=float(weight*dx**3))
