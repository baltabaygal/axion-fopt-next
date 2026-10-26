"""Stored energy around a translating pocket; not a moving-boundary flux.

Uses symmetric bond energies, a periodic unwrapped false-vacuum centroid,
and the same radius as the fixed monitor. This diagnoses whether translation
has carried bound material across that fixed sphere. No lifetime or rest-mass
claim follows from a few snapshots.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import dilaton_potential,radius,read_state
from pocket_geometry import periodic_centroid


def moving_budget(state,box,monitor_radius,lam,dv,mass2,power,vacuum_stiffening=0.):
    n=len(state);dx=box/n;phi,a,pphi,pa=(state[...,j] for j in range(4));fv=phi<.5
    centre,coherence=periodic_centroid(fv,box)
    distance=radius(n,box,centre);inside=distance<monitor_radius
    phi_density=.5*pphi**2+dilaton_potential(phi,lam,dv,vacuum_stiffening)
    a_density=.5*pa**2+.5*mass2*phi**power*a*a
    momentum=[];global_momentum=[]
    for axis in range(3):
        gradient_phi=(np.roll(phi,-1,axis)-np.roll(phi,1,axis))/(2*dx)
        gradient_a=(np.roll(a,-1,axis)-np.roll(a,1,axis))/(2*dx)
        momentum_density=-pphi*gradient_phi-pa*gradient_a
        momentum.append(float(np.sum(momentum_density[inside])*dx**3))
        global_momentum.append(float(np.sum(momentum_density)*dx**3))
        phi_density+=((np.roll(phi,-1,axis)-phi)**2+(phi-np.roll(phi,1,axis))**2)/(4*dx*dx)
        a_density+=((np.roll(a,-1,axis)-a)**2+(a-np.roll(a,1,axis))**2)/(4*dx*dx)
    e_phi=float(np.sum(phi_density[inside])*dx**3);e_a=float(np.sum(a_density[inside])*dx**3)
    return dict(centroid=centre.tolist(),centroid_coherence=coherence,
        maximum_fv_radius_about_centroid=float(distance[fv].max()),
        minimum_monitor_clearance=float(monitor_radius-distance[fv].max()),
        stored_phi_energy=e_phi,stored_axion_energy=e_a,stored_total_energy=e_phi+e_a,
        global_total_energy=float(np.sum(phi_density+a_density)*dx**3),
        stored_momentum=momentum,global_momentum=global_momentum,
        interpretation='Instantaneous storage, not cumulative radiation or a moving-volume continuity result.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((args.case/'input.json').read_text())
    if cfg.get('expansion'):
        raise ValueError('This audit supports static fields only')
    paths=[args.case/'native/initial.bin']+sorted((args.case/'native').glob('frame_*.bin'))+[args.case/'native/final.bin']
    result={}
    for path in paths:
        meta,s=read_state(path)
        if str(meta['u']) in result:
            continue
        result[str(meta['u'])]=dict(time=meta['u'],snapshot=str(path.resolve()),
            source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            **moving_budget(s,meta['box'],cfg['monitor_radius'],cfg['lambda'],cfg['delta_v'],cfg['mass_tv']**2,cfg.get('mass_power',2),cfg.get('vacuum_stiffening',0.)))
        del s
    parent=Path(__file__).parent
    report=dict(snapshots=result,sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest()
                    for f in ('moving_snapshot_budget.py','pocket_geometry.py','model.py')})
    (args.output/'budget.json').write_text(json.dumps(report,indent=2)+'\n')
    for t,r in result.items():
        print(t,'centroid',r['centroid'],'stored total',r['stored_total_energy'],'clearance',r['minimum_monitor_clearance'])


if __name__=='__main__':
    main()
