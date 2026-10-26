"""Slab-wise static snapshot accounting with bounded temporary arrays.

Preserves the legacy forward-bond sector allocation and arithmetic FV
centroid for comparison. That centroid is unsuitable for a pocket crossing
the periodic cut; use the separate moving-volume monitor for retention.
Memory mapping avoids a full Python snapshot copy, not OS page caching.
"""
import hashlib
import math
from pathlib import Path
import struct

import numpy as np
from model import MAX_N,dilaton_potential


def budget(path,lam=2.,dv=.01,mass=2.,margin=6.,mass_power=2,vacuum_stiffening=0.,slab=8):
    path=Path(path)
    if not isinstance(slab,int) or slab<1 or mass_power not in (2,4):
        raise ValueError('Require positive integer slab and mass power 2/4')
    with path.open('rb') as f:
        header=f.read(28)
        if len(header)!=28 or header[:8]!=b'PKST0001':
            raise ValueError('Invalid snapshot header')
        n,box,time=struct.unpack('<I2d',header[8:])
        if not 4<=n<=MAX_N or not np.isfinite([box,time]).all() or box<=0 or path.stat().st_size!=28+32*n**3:
            raise ValueError('Invalid snapshot metadata/size')
        f.seek(0);digest=hashlib.file_digest(f,'sha256').hexdigest()
    s=np.memmap(path,mode='r',dtype='<f8',offset=28,shape=(n,n,n,4))
    dx=box/n;cell=dx**3;x=(np.arange(n)+.5)*dx
    count=0;moments=np.zeros(3)
    for start in range(0,n,slab):
        end=min(n,start+slab);block=s[start:end]
        if not np.isfinite(block).all():raise ValueError('Nonfinite snapshot')
        fv=block[...,0]<.5;count+=int(fv.sum())
        moments+=np.array([np.sum(fv.sum(axis=(1,2))*x[start:end]),
                           np.sum(fv.sum(axis=(0,2))*x),np.sum(fv.sum(axis=(0,1))*x)])
    radius=np.cbrt(3*count*cell/(4*np.pi));centre=moments/count if count else None
    sums={key:[] for key in ('dilaton_kinetic','dilaton_gradient','dilaton_potential',
          'axion_kinetic','axion_gradient','coupling','axion_fv','axion_wall','axion_tv','axion_far')}
    far_phi=[];exchanges=[]
    for start in range(0,n,slab):
        end=min(n,start+slab);block=s[start:end]
        phi,a,pphi,pa=(block[...,j] for j in range(4))
        gp=np.zeros_like(phi);ga=np.zeros_like(a)
        for axis in range(3):
            following=s[(np.arange(start,end)+1)%n] if axis==0 else np.roll(block,-1,axis)
            gp+=(following[...,0]-phi)**2/(2*dx*dx)
            ga+=(following[...,1]-a)**2/(2*dx*dx)
        kp=.5*pphi**2;ka=.5*pa**2
        up=dilaton_potential(phi,lam,dv,vacuum_stiffening)
        coupling=.5*mass**2*phi**mass_power*a*a
        axion=ka+ga+coupling;fv=phi<.5;tv=phi>.99
        values=dict(dilaton_kinetic=kp,dilaton_gradient=gp,dilaton_potential=up,
                    axion_kinetic=ka,axion_gradient=ga,coupling=coupling,
                    axion_fv=axion[fv],axion_wall=axion[~fv & ~tv],axion_tv=axion[tv])
        far=np.zeros_like(fv)
        if count:
            ds=[np.abs(v-c) for v,c in zip((x[start:end],x,x),centre)]
            ds=[np.minimum(v,box-v) for v in ds]
            far=np.sqrt(ds[0][:,None,None]**2+ds[1][None,:,None]**2+ds[2][None,None,:]**2)>radius+margin
            far_phi.append(float(np.sum((kp+gp+up)[far]))*cell)
        values['axion_far']=axion[far]
        for key,value in values.items():sums[key].append(float(np.sum(value))*cell)
        exchanges.append(float(np.sum(.5*mass_power*mass**2*phi**(mass_power-1)*a*a*pphi))*cell)
    e={key:math.fsum(value) for key,value in sums.items()}
    e.update(fv_volume=count*cell,fv_radius=float(radius))
    e['total']=math.fsum(e[key] for key in ('dilaton_kinetic','dilaton_gradient','dilaton_potential','axion_kinetic','axion_gradient','coupling'))
    return dict(time=time,energies=e,source_sha256=digest,mass_power=mass_power,
                instantaneous_axion_exchange=math.fsum(exchanges),
                far_dilaton_energy=math.fsum(far_phi) if count else None,
                far_boundary_radius=float(radius+margin) if count else None)
