"""Stationary second-wall displacement with a fixed local incident packet."""
import numpy as np
from planar_initial import PlanarSeed,modes
from paper_reference import Parameters


def extended_seed(dx=.125,speed=0.,amplitude=.001):
    if not np.isfinite(dx) or dx<=0 or not np.isfinite(speed) or not 0<=speed<1:
        raise ValueError("Invalid grid or wall speed")
    box=768.;origin=-256.;n=round(box/dx)
    gamma=1/np.sqrt(1-speed*speed)
    if not np.isfinite(dx) or dx<=0 or abs(n*dx-box)>1e-10 or not 8<=n<=16384:
        raise ValueError('Invalid extended planar grid')
    physical=(np.arange(n)+.5)*dx+origin
    s=np.zeros((n,4));s[:,0]=1.
    for image in (-box,0.,box):
        for location,sign,velocity in ((256.,1.,-speed),(origin,-1.,speed)):
            t=np.tanh(gamma*(physical-location-image))
            s[:,0]+=.5*sign*t
            s[:,2]+=-velocity*.5*sign*gamma*(1-t*t)
    seed=PlanarSeed(amplitude=amplitude);p=Parameters(dv=0.);indices,k,omega,weights=modes(seed,p)
    for kk,ww,weight in zip(k,omega,weights):
        phase=kk*(physical-seed.centre)
        s[:,1]+=seed.amplitude*weight*np.cos(phase)
        s[:,3]+=seed.amplitude*weight*ww*np.sin(phase)
    untapered=s[:,[1,3]].copy()
    # Keep all coefficients and fields unchanged in the common physical core.
    # Suppress only far tails so the enlarged domain is exactly periodic.
    z=np.clip((abs(physical-seed.centre)-176.)/48.,0.,1.)
    cutoff=1-z**3*(10+z*(-15+6*z))
    s[:,1]*=cutoff;s[:,3]*=cutoff
    return physical,s,untapered


def moving_final_core(speed,duration,width=1.):
    """Interior excluding eight widths from the known nominal wall paths.

    Initial walls are at 0 and 256 in the original physical coordinates.
    This is a comparison region, not a measured causal or material boundary.
    """
    if not np.isfinite([speed,duration,width]).all() or not 0<=speed<1 or duration<0 or width<=0:
        raise ValueError('Invalid moving comparison geometry')
    lower=speed*duration+8*width;upper=256-speed*duration-8*width
    if lower>=upper:raise ValueError('No common interior remains between nominal wall paths')
    return lower,upper
