"""Deterministic continuum-defined initial data for planar coupled fields.

Only O(N) arrays are allocated. Preparing these data does not run a solver.
Both wall velocities and the positive-k packet momentum are explicit.
"""
from dataclasses import dataclass,asdict
import math
import numpy as np
from paper_reference import Parameters,terms,head_on_reflection


@dataclass(frozen=True)
class PlanarSeed:
    box:float=512.
    width:float=1.
    speed:float=0.
    centre:float=128.
    envelope_width:float=16.
    omega:float=.5
    amplitude:float=1e-3
    sigma_cut:float=6.

    def validate(self,p):
        if (not np.isfinite(list(asdict(self).values())).all() or self.box<=0 or self.width<=0 or
            not 0<=self.speed<1 or not 0<self.centre<self.box/2 or self.envelope_width<=0 or
            self.omega<=p.mass_false or self.amplitude<0 or self.sigma_cut<=0):
            raise ValueError('Invalid planar seed')
        if p.dv!=0. or abs(p.lam-8./self.width**2)>1e-12*p.lam:
            raise ValueError('This wall seed requires degenerate quartic and lambda=8/width^2')


def modes(seed,p):
    seed.validate(p);k0=math.sqrt(seed.omega**2-p.mass_false**2)
    low=max(1,math.ceil((k0-seed.sigma_cut/seed.envelope_width)*seed.box/(2*math.pi)))
    high=math.floor((k0+seed.sigma_cut/seed.envelope_width)*seed.box/(2*math.pi))
    if high<low:raise ValueError('No packet modes in this box')
    indices=np.arange(low,high+1);k=2*math.pi*indices/seed.box
    weights=np.exp(-.5*(seed.envelope_width*(k-k0))**2);weights/=weights.sum()
    return indices,k,np.sqrt(k*k+p.mass_false**2),weights


def prepare(n,seed=PlanarSeed(),p=Parameters(dv=0.)):
    if not isinstance(n,int) or n<8 or n>16384:raise ValueError('Require N=8..16384')
    indices,k,omega,weights=modes(seed,p)
    if indices[-1]>=n//2:raise ValueError('Packet crosses Nyquist; increase N')
    dx=seed.box/n;x=(np.arange(n)+.5)*dx;gamma=1/math.sqrt(1-seed.speed**2)
    state=np.zeros((n,4));state[:,0]=1.
    for image in (-seed.box,0.,seed.box):
        for location,sign,velocity in ((seed.box/2,1.,-seed.speed),(0.,-1.,seed.speed)):
            z=gamma*(x-location-image)/seed.width;t=np.tanh(z)
            state[:,0]+=.5*sign*t
            state[:,2]+=-velocity*.5*sign*gamma/seed.width*(1-t*t)
    # Identical continuum coefficients on every grid. In particular, do not
    # retune pi_a to lattice dispersion when comparing spatial resolutions.
    for kk,ww,weight in zip(k,omega,weights):
        phase=kk*(x-seed.centre)
        state[:,1]+=seed.amplitude*weight*np.cos(phase)
        state[:,3]+=seed.amplitude*weight*ww*np.sin(phase)
    return x,state


def audit(n,seed=PlanarSeed(),p=Parameters(dv=0.)):
    x,s=prepare(n,seed,p);dx=seed.box/n;indices,k,omega,weights=modes(seed,p)
    va=terms(s[:,0],s[:,1],p)[1]
    gradient=(np.roll(s[:,1],-1)-s[:,1])/dx
    density=.5*s[:,3]**2+va+.25*(gradient**2+np.roll(gradient,1)**2)
    total=float(density.sum()*dx)
    ak=np.fft.fft(s[:,1])/n;pk=np.fft.fft(s[:,3])/n
    def left_fraction(frequencies):
        right=.5*(ak[indices]+1j*pk[indices]/frequencies)
        left=.5*(ak[indices]-1j*pk[indices]/frequencies)
        er=np.sum(frequencies**2*abs(right)**2);el=np.sum(frequencies**2*abs(left)**2)
        return float(el/(er+el)) if er+el else 0.
    lattice_omega=np.sqrt(p.mass_false**2+(2*np.sin(k*dx/2)/dx)**2)
    modal_energy=omega**2*weights**2
    wall_omega=np.array([head_on_reflection(float(ww),p.mass_false,seed.speed)['wall_frame_frequency'] for ww in omega])
    mt=p.mass_false*np.exp(p.kappa/2)
    shifted=lambda xx:(xx+seed.box/2)%seed.box-seed.box/2
    distances=np.minimum(abs(shifted(x)),abs(shifted(x-seed.box/2)))
    phi=s[:,0]
    lap=(np.roll(phi,-1)+np.roll(phi,1)-2*phi)/dx**2
    vacuum_force=p.lam/2*phi*(1-phi)*(1-2*phi)
    residual=(1-seed.speed**2)*lap-vacuum_force
    return dict(n=n,dx=dx,seed=asdict(seed),potential_parameters=asdict(p),modes=len(k),
        k_band=[float(k.min()),float(k.max())],omega_band=[float(omega.min()),float(omega.max())],
        initial_axion_energy_per_area=total,
        continuum_left_moving_fraction=left_fraction(omega),lattice_left_moving_fraction=left_fraction(lattice_omega),
        initial_energy_in_tv_fraction=float(density[phi>.9].sum()*dx/total) if total else 0.,
        initial_energy_near_walls_fraction=float(density[distances<8*seed.width].sum()*dx/total) if total else 0.,
        initial_modal_energy_above_wall_threshold=float(modal_energy[wall_omega>=mt].sum()/modal_energy.sum()),
        boosted_wall_10_90_cells=float(2*np.arctanh(.8)*seed.width*np.sqrt(1-seed.speed**2)/dx),
        continuum_kink_discrete_residual=float(np.max(np.abs(residual))),
        centre_arrival_estimate=(seed.box/2-seed.centre)/(math.sqrt(seed.omega**2-p.mass_false**2)/seed.omega+seed.speed),
        other_wall_signal_to_target_estimate=(seed.box/2)/(1+seed.speed),
        interpretation='Continuum seed sampled on a grid. Arrival/image estimates ignore noncompact tails and subsequent acceleration; empty-wall and boundary controls required.')
