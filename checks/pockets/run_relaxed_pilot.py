"""Matched-action, field-only preparation pilot; see RELAXED_PILOT.md."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from model import EXECUTABLE, config, grid, history, read_slice, read_state, real_sph_harm, run
from relax_reference import RadialEnergy
from wall_diagnostics import local_widths


def field_state(n, box, rgrid, phi, modes, seed=0):
    x=grid(n,box)-box/2
    X,Y,Z=(v.ravel() for v in np.meshgrid(x,x,x,indexing='ij'))
    r=np.sqrt(X*X+Y*Y+Z*Z)
    state=np.zeros((n**3,4))
    state[:,0]=np.interp(r,rgrid,phi,left=phi[0],right=1.)
    keep=r<30.; rk=r[keep]
    theta=np.arccos(np.clip(Z[keep]/rk,-1.,1.));azimuth=np.arctan2(Y[keep],X[keep])
    a=np.zeros(len(rk));pa=np.zeros(len(rk));rng=np.random.default_rng(seed)
    for l in sorted({m['l'] for m in modes}):
        harmonics=[real_sph_harm(l,m,theta,azimuth) for m in range(-l,l+1)]
        for mode in (m for m in modes if m['l']==l):
            u=mode['u'].copy()
            if u[np.argmax(np.abs(u))]<0:
                u=-u
            radial=np.interp(rk,rgrid,u/rgrid,right=0.)
            for y in harmonics:
                phase=rng.uniform(0,2*np.pi); action=mode['action']; omega=mode['omega']
                a+=np.sqrt(2*action/omega)*np.cos(phase)*radial*y
                pa+=np.sqrt(2*action*omega)*np.sin(phase)*radial*y
    state[keep,1]=a;state[keep,3]=pa
    return state.reshape(n,n,n,4)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,default=Path('outputs/pockets/relaxed-reference-h005-01'))
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    report=json.loads((args.reference/'report.json').read_text())
    if not report['reference_pass']:
        raise ValueError('Preparation reference must pass')
    reference=json.loads((args.reference/'contract.json').read_text())
    groups={int(l):m for l,m in reference['groups'].items()}
    system=RadialEnergy(reference['h'],reference['rmax'],groups)
    data=np.load(args.reference/'profile.npz')
    profiles=dict(tanh=data['initial_phi'],relaxed=data['phi'])
    r=data['r'];modes={name:system.evaluate(phi,return_modes=True)[2]['modes'] for name,phi in profiles.items()}
    args.output.mkdir(parents=True,exist_ok=False)
    parent=Path(__file__).parent
    manifest=dict(reference=str(args.reference.resolve()),seed=0,
                  states={name:dict(continuum_gas_energy=sum((2*m['l']+1)*m['action']*m['omega'] for m in group),
                                    highest_omega=max(m['omega'] for m in group)) for name,group in modes.items()},
                  sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest() for f in (
                      'RELAXED_PILOT.md','run_relaxed_pilot.py','relax_reference.py','pocket_model.h','evolve.cpp')},
                  reference_profile_sha256=hashlib.sha256((args.reference/'profile.npz').read_bytes()).hexdigest(),
                  executable_sha256=hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest())
    (args.output/'contract.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def execute(name):
        s=field_state(128,64.,r,profiles[name],modes[name])
        cfg=config(128,64.,.0125,80.,lam=2.,dv=.01,mass_tv=2.,every=40,slices=320)
        error=None
        try:
            run(s,cfg,args.output/name)
        except RuntimeError as exc:
            error=str(exc)
        h=history(args.output/name)
        _,transferred=read_state(args.output/name/'native/initial.bin')
        e=sum(h[k] for k in ('axion_kinetic','axion_gradient','coupling'))
        drift=float(np.max(np.abs(h['total']/h['total'][0]-1)))
        walls=[]
        for path in sorted((args.output/name/'native').glob('slice_*.bin')):
            meta,plane=read_slice(path)
            walls.append(dict(time=meta['u'],**local_widths(plane[...,0],.5)))
        gates=dict(transfer=bool(np.array_equal(s,transferred)),completed=error is None and bool(np.isclose(h['u'][-1],80.)),
                   energy=drift<5e-4)
        outcome=dict(gates=gates,energy_drift=drift,solver_error=error,
                     initial_radius=float(h['fv_radius'][0]),final_radius=float(h['fv_radius'][-1]),
                     radius_loss=float(h['fv_radius'][0]-h['fv_radius'][-1]),
                     initial_axion_energy=float(e[0]),axion_depletion=float(1-e[-1]/e[0]),
                     final_far_storage=float(h['axion_far'][-1]/e[0]),
                     walls=walls,slice_diagnostic_pass=all(w['diagnostic_pass'] for w in walls))
        (args.output/name/'summary.json').write_text(json.dumps(outcome,indent=2,allow_nan=False)+'\n')
        print(name,{k:v for k,v in outcome.items() if k!='walls'},flush=True)
        return outcome
    with ThreadPoolExecutor(2) as pool:
        outcomes=dict(zip(profiles,pool.map(execute,profiles)))
    a,b=outcomes['tanh'],outcomes['relaxed']
    improved=bool(a['radius_loss']>0 and a['axion_depletion']>0
                  and b['radius_loss']<=.5*a['radius_loss'] and b['axion_depletion']<=.5*a['axion_depletion']
                  and all(all(o['gates'].values()) and o['slice_diagnostic_pass'] for o in outcomes.values()))
    result=dict(runs=outcomes,improved=improved,interpretation='Field-only pilot, not long survival or continuum proof.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
