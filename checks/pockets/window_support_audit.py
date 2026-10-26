"""Slab-wise energy and threshold geometry near a moving window edge.

Threshold sets include any disconnected waves. This does not identify bound
energy or irreversible radiation; it checks sensitivity to window support.
"""
import math
import numpy as np
from model import dilaton_potential


def audit(state,box,centre,inner,outer,lam,dv,mass2,power=4,b=0.,slab=8):
    n=len(state);dx=box/n;cell=dx**3
    if state.shape!=(n,n,n,4) or not 0<inner<outer<box/2 or not isinstance(slab,int) or slab<1:
        raise ValueError('Invalid field/window/slab')
    x=(np.arange(n)+.5)*dx
    ds=[(x-c+box/2)%box-box/2 for c in centre]
    sums={k:[] for k in ('weighted_phi','weighted_axion','inner_phi','inner_axion',
                         'annulus_phi','annulus_axion','outside_phi','outside_axion',
                         'interface_outside_inner_phi','interface_outside_inner_axion')}
    extents={str(v):None for v in (.1,.5,.9,.99)}
    for start in range(0,n,slab):
        end=min(n,start+slab);block=state[start:end]
        if not np.isfinite(block).all():raise ValueError('Nonfinite fields')
        phi,a,pphi,pa=(block[...,j] for j in range(4))
        ep=.5*pphi**2+dilaton_potential(phi,lam,dv,b)
        ea=.5*pa**2+.5*mass2*phi**power*a*a
        for axis in range(3):
            plus=state[(np.arange(start,end)+1)%n] if axis==0 else np.roll(block,-1,axis)
            minus=state[(np.arange(start,end)-1)%n] if axis==0 else np.roll(block,1,axis)
            ep+=((plus[...,0]-phi)**2+(phi-minus[...,0])**2)/(4*dx*dx)
            ea+=((plus[...,1]-a)**2+(a-minus[...,1])**2)/(4*dx*dx)
        r=np.sqrt(ds[0][start:end,None,None]**2+ds[1][None,:,None]**2+ds[2][None,None,:]**2)
        s=np.clip((r-inner)/(outer-inner),0.,1.);weight=1-s**3*(10+s*(-15+6*s))
        for level in extents:
            selected=r[phi<float(level)]
            if len(selected):extents[level]=max(extents[level] or 0.,float(selected.max()))
        masks=dict(inner=r<inner,annulus=(r>=inner)&(r<outer),outside=r>=outer,
                   interface_outside_inner=(r>=inner)&(phi>.1)&(phi<.9))
        for sector,density in (('phi',ep),('axion',ea)):
            sums['weighted_'+sector].append(float(np.sum(weight*density))*cell)
            for region,mask in masks.items():sums[region+'_'+sector].append(float(np.sum(density[mask]))*cell)
    result={k:math.fsum(v) for k,v in sums.items()}
    for region in ('weighted','inner','annulus','outside','interface_outside_inner'):
        result[region+'_total']=result[region+'_phi']+result[region+'_axion']
    result['threshold_max_radii']=extents
    return result
