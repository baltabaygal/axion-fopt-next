"""Static snapshot energy-channel audit using the compiled model's stencil.

Stored far energies are not cumulative radiation. The instantaneous exchange
Q=integral m_TV^2 phi a^2 phi_dot is signed positive into the axion sector.
This allocation assigns the entire coupling potential to the axion sector.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import dilaton_potential, energies, radius, read_state


def budget(path, lam=2., dv=.01, mass=2., margin=6., mass_power=2, vacuum_stiffening=0.):
    meta,s=read_state(path);n=meta['n'];box=meta['box'];dx=box/n
    phi,a,pphi,pa=(s[...,j] for j in range(4))
    e=energies(s,box,lam,dv,mass**2,far_margin=margin,mass_power=mass_power,vacuum_stiffening=vacuum_stiffening)
    fv=phi<.5
    result=dict(time=meta['u'],energies=e,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                mass_power=mass_power,
                instantaneous_axion_exchange=float(np.sum(.5*mass_power*mass**2*phi**(mass_power-1)*a*a*pphi)*dx**3))
    if fv.any():
        x=(np.arange(n)+.5)*dx
        centre=[float(x[i].mean()) for i in np.nonzero(fv)]
        far=radius(n,box,centre)>e['fv_radius']+margin
        density=.5*pphi**2+dilaton_potential(phi,lam,dv,vacuum_stiffening)
        for axis in range(3):
            density+=(np.roll(phi,-1,axis)-phi)**2/(2*dx**2)
        result['far_dilaton_energy']=float(np.sum(density[far])*dx**3)
        result['far_boundary_radius']=e['fv_radius']+margin
    else:
        result['far_dilaton_energy']=None
        result['far_boundary_radius']=None
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshots',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mass-power',type=int,choices=(2,4),default=2)
    args=p.parse_args()
    if args.output.exists():
        p.error('Output must be fresh')
    result=dict(snapshots={str(path.resolve()):budget(path,mass_power=args.mass_power) for path in args.snapshots},
                method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                interpretation='Static phi^p model snapshot energy allocation; far storage is not irreversible loss.')
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'budget.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for path,row in result['snapshots'].items():
        print(path,'far dilaton',row['far_dilaton_energy'],'exchange',row['instantaneous_axion_exchange'])


if __name__=='__main__':
    main()
