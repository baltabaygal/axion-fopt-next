"""Tiny compiled spatial-order test using a controlled weak free mode."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from paper_reference import Parameters
from planar_stencils import symbol

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    exe=args.executable.resolve();build=json.loads((exe.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale spatial build')
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Binary mismatch')
    if build.get('supports_line_spatial_orders')!=[2,4]:raise RuntimeError('Both spatial orders required')
    gates=dict(second_min_ratio=3.,fourth_min_ratio=12.,semidiscrete_scaled_error=1e-8,maximum_solver_seconds=5.)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    p=Parameters();box=8.;dt=.00025;duration=.1;amplitude=1e-6;k=2*np.pi/box
    omega=np.sqrt(k*k+p.mass_false**2);rows=[];spent=0.
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    for order in (2,4):
        for n in (16,32,64):
            x=(np.arange(n)+.5)*box/n;s=np.zeros((n,4));s[:,1]=amplitude*np.cos(k*x);s[:,3]=amplitude*omega*np.sin(k*x)
            directory=args.output/f'order{order}-n{n}';directory.mkdir();np.savetxt(directory/'input.txt',s,fmt='%.17g')
            params=dict(N=n,lSide=box,dt=dt,tMax=duration,expansion=False,evolver='PV2',PS_type=2,line_spatial_order=order,
                **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,output_every=400,
                initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
            (directory/'input.in').write_text(''.join(f'{key} = {str(v).lower() if isinstance(v,bool) else v}\n' for key,v in params.items()))
            start=time.monotonic();result=subprocess.run([str(exe),'input='+str((directory/'input.in').resolve())],env=env,capture_output=True,text=True,timeout=max(.1,gates['maximum_solver_seconds']-spent))
            spent+=time.monotonic()-start;(directory/'run.log').write_text(result.stdout+result.stderr)
            if result.returncode:raise RuntimeError(result.stderr[-1000:])
            final=np.loadtxt(directory/'native/final.txt');wl=np.sqrt(p.mass_false**2+symbol(k,box/n,order))
            reference=np.column_stack((s[:,1]*np.cos(wl*duration)+s[:,3]/wl*np.sin(wl*duration),
                -s[:,1]*wl*np.sin(wl*duration)+s[:,3]*np.cos(wl*duration)))
            continuum=np.column_stack((amplitude*np.cos(k*x-omega*duration),amplitude*omega*np.sin(k*x-omega*duration)))
            scales=np.array([amplitude,amplitude*omega]);np.savetxt(directory/'semidiscrete-reference.txt',reference,fmt='%.17g')
            rows.append(dict(order=order,n=n,semidiscrete_scaled_error=float(np.max(abs(final[:,[1,3]]-reference)/scales)),
                             continuum_scaled_error=float(np.max(abs(final[:,[1,3]]-continuum)/scales)),
                             induced_phi_max=float(np.max(abs(final[:,0])))))
    ratios={str(order):[r[i]['continuum_scaled_error']/r[i+1]['continuum_scaled_error'] for i in (0,1)]
            for order in (2,4) for r in [[v for v in rows if v['order']==order]]}
    checks=dict(semidiscrete=all(r['semidiscrete_scaled_error']<gates['semidiscrete_scaled_error'] for r in rows),
        second_order=all(v>gates['second_min_ratio'] for v in ratios['2']),fourth_order=all(v>gates['fourth_min_ratio'] for v in ratios['4']))
    report=dict(kind=__doc__,scope='Canonical periodic line, weak-mode linear reference; not pocket physics',
        box=box,dt=dt,duration=duration,time='Minkowski cosmic time',amplitude=amplitude,k=k,omega=omega,
        nonlinear_assumption='Induced dilaton is measured; coupling correction to axion mode is O(amplitude^2)',
        gates=gates,rows=rows,ratios=ratios,checks=checks,passed=all(checks.values()),total_solver_seconds=spent,build=build,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('validate_planar_spatial_order.py','planar_stencils.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({key:report[key] for key in ('rows','ratios','checks','passed','total_solver_seconds')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
