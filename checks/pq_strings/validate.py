"""Predeclared correctness gates for the compiled two-real-component PQ model."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from model import (ROOT, EXECUTABLE, config, run, read_state, reference, potential,
                   force, deterministic, defects, energies, coefficients)

GATES = dict(force_difference=2e-7, vacuum_mass_error=2e-6, final_ode_error=2e-4,
             refinement_ratio=3.5, transfer_error=0., energy_snapshot_error=2e-12,
             energy_balance_relative=2e-4)

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=EXECUTABLE);args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'contract.json').write_text(json.dumps(GATES,indent=2)+'\n')
    checks={};details={}
    rng=np.random.default_rng(0);x=rng.normal(size=100);y=rng.normal(size=100);eps=1e-6
    numerical=np.stack(((potential(x+eps,y,.4,.03)-potential(x-eps,y,.4,.03))/(2*eps),
                        (potential(x,y+eps,.4,.03)-potential(x,y-eps,.4,.03))/(2*eps)),axis=-1)
    err=float(np.max(np.abs(numerical-force(x,y,.4,.03))))
    details['force_difference']=err;checks['force_difference']=err<GATES['force_difference']
    checks['origin_limit']=bool(np.array_equal(force(np.array(0.),np.array(0.),.4,.03),[0.,0.]))
    hessian=np.stack(((force(np.array(1+eps),np.array(0.),.4,.03)-force(np.array(1-eps),np.array(0.),.4,.03))/(2*eps),
                      (force(np.array(1.),np.array(eps),.4,.03)-force(np.array(1.),np.array(-eps),.4,.03))/(2*eps)))
    mass_error=float(np.max(np.abs(hessian-np.diag([.8,.03]))));details['vacuum_mass_error']=mass_error
    checks['vacuum_masses']=mass_error<GATES['vacuum_mass_error']
    d=deterministic(16,8.);topology,w=defects(d,8.,10.)
    checks['known_winding']=bool(topology['plaquettes']==4*16 and topology['signed_winding']==[0,0,0]
                                and np.count_nonzero(w[0][...,0])==4)
    details['known_winding']=topology
    for case in ('homogeneous_static','inhomogeneous_static','inhomogeneous_radiation'):
        n=4;box=4.;s=np.zeros((n,n,n,4));s[...,0]=1.05;s[...,1]=.15;s[...,2]=.02;s[...,3]=-.03
        expansion=case.endswith('radiation')
        if expansion:
            grid=(np.arange(n)+.5)*box/n
            X,Y,Z=np.meshgrid(grid,grid,grid,indexing='ij')
            s[...,0]+=.1*np.sin(2*np.pi*X/box);s[...,1]+=.07*np.cos(2*np.pi*Y/box)
            s[0,0,0,:2]=0. # Exercise the finite Cartesian force at an exact core.
        if case == 'inhomogeneous_static':
            grid=(np.arange(n)+.5)*box/n
            X,Y,Z=np.meshgrid(grid,grid,grid,indexing='ij')
            s[...,0]=.8*np.cos(2*np.pi*X/box)+.1*np.sin(2*np.pi*Z/box)
            s[...,1]=.8*np.sin(2*np.pi*X/box)+.07*np.cos(2*np.pi*Y/box)
            s[0,0,0,:2]=0.
        errors=[]
        for dt in (.005,.0025):
            cfg=config(n,box,dt,1.,lam=.4,origin=2.,star=2.5,chi_star=.03,
                       exponent=4.,expansion=expansion,every=1)
            directory=args.output/f'{case}_{dt}'
            cfg=run(s,cfg,directory,args.executable)
            h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True)
            ref=reference(s,cfg,[1.])[0];meta,final=read_state(directory/'native/final.bin')
            error=float(np.max(np.abs(final-ref)));errors.append(error)
            _,initial=read_state(directory/'native/initial.bin')
            checks[f'{case}_{dt}_transfer']=bool(np.array_equal(initial,s))
            tau,a,chi=coefficients(cfg,1.)
            checks[f'{case}_{dt}_time']=abs(meta['tau']-tau)<1e-14 and abs(meta['a']-a)<1e-14 and abs(meta['chi']-chi)<1e-14
            independent=energies(final,box,a,.4,chi)
            energy_err=max(abs(independent[k]-h[k][-1]) for k in independent)
            checks[f'{case}_{dt}_energy']=energy_err<GATES['energy_snapshot_error']
            work=np.trapezoid(h['source']+h['dilution'],h['u'])
            balance=abs(h['total'][-1]-h['total'][0]-work)/max(abs(h['total'][0]),abs(h['total'][-1]),1e-12)
            checks[f'{case}_{dt}_balance']=balance<GATES['energy_balance_relative']
            details[f'{case}_{dt}']=dict(final_ode_error=error,energy_error=energy_err,balance_relative=float(balance))
        checks[case+'_ode']=max(errors)<GATES['final_ode_error']
        checks[case+'_refinement']=errors[0]/errors[1]>GATES['refinement_ratio']
        details[case+'_ratio']=errors[0]/errors[1]
    # Bad state and unsupported configuration must fail without creating native output.
    import subprocess
    bad=args.output/'bad.in'
    base=(args.output/'homogeneous_static_0.005/input.in').read_text()
    for label,replacement in [('nonfinite','lambda = nan'),('fractional_grid','N = 4.5'),('wrong_background','H0 = 2')]:
        text=base.replace('expansion = false','expansion = true')
        key=replacement.split(' = ')[0]
        text='\n'.join(replacement if line.startswith(key+' = ') else line for line in text.splitlines())
        text=text.replace(str((args.output/'homogeneous_static_0.005/native').resolve()),str((args.output/label).resolve()))
        bad.write_text(text)
        result=subprocess.run([str(args.executable),'input='+str(bad.resolve())],capture_output=True,text=True)
        checks['reject_'+label]=result.returncode!=0 and not (args.output/label).exists()
    checks = {k: bool(v) for k, v in checks.items()}
    report=dict(gates=GATES,checks=checks,details=details,passed=all(checks.values()),
                executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest())
    (args.output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed'] else 1)
if __name__=='__main__':main()
