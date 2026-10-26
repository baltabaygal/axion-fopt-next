"""Short one-thread wall controls, with a strict total runtime budget."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from planar_diagnostics import wall_geometry

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--prepared',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale line executable')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    gates=dict(duration=4.,maximum_total_solver_seconds=12.,position_drift=.02,phase_velocity_change=.005,
               width_relative_change=.05,relative_energy_drift=2e-4,finest_continuum_field_error=.002)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    rows=[];solver_seconds=0.
    for speed in (0.,.3):
        for n in (2048,4096):
            source=args.prepared/f'u{speed}-n{n}-empty-wall';params=json.loads((source/'parameters.json').read_text())
            directory=args.output/f'u{speed}-n{n}';directory.mkdir()
            initial=np.loadtxt(source/'input.txt');np.savetxt(directory/'input.txt',initial,fmt='%.17g')
            params.update(tMax=4.,initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
            (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
            (directory/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
            remaining=gates['maximum_total_solver_seconds']-solver_seconds
            if remaining<=0:raise RuntimeError('Short-control runtime budget exhausted; no further runs launched')
            start=time.monotonic()
            run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],env=environment,
                               capture_output=True,text=True,timeout=remaining)
            elapsed=time.monotonic()-start;solver_seconds+=elapsed
            (directory/'run.log').write_text(run.stdout+run.stderr)
            if run.returncode:raise RuntimeError(run.stderr[-2000:])
            final=np.loadtxt(directory/'native/final.txt');h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True)
            box=params['lSide'];before=wall_geometry(initial,box);after=wall_geometry(final,box)
            expected=(before['position']-speed*4.)%box
            displacement=(after['position']-expected+box/2)%box-box/2
            x=(np.arange(n)+.5)*box/n;gamma=1/np.sqrt(1-speed**2);ref=np.zeros((n,4));ref[:,0]=1.
            for image in (-box,0.,box):
                for location,sign,velocity in ((box/2,1.,-speed),(0.,-1.,speed)):
                    z=gamma*(x-location-velocity*4.-image);t=np.tanh(z)
                    ref[:,0]+=.5*sign*t;ref[:,2]+=-velocity*.5*sign*gamma*(1-t*t)
            np.savetxt(directory/'continuum-final.txt',ref,fmt='%.17g')
            row=dict(speed=speed,n=n,elapsed_seconds=elapsed,before=before,after=after,
                position_drift_from_continuum=float(displacement),phase_velocity_change=abs(after['phase_velocity']-before['phase_velocity']),
                width_relative_change=abs(after['width_10_90']/before['width_10_90']-1.),
                relative_energy_drift=float(np.max(np.abs(h['total']-h['total'][0]))/h['total'][0]),
                continuum_field_error=float(np.max(np.abs(final-ref))),
                zero_axion_exact=bool(np.all(final[:,[1,3]]==0.) and np.all(h['axion_energy']==0.) and np.all(h['exchange']==0.)))
            row['passed']=bool(abs(row['position_drift_from_continuum'])<gates['position_drift'] and
                row['phase_velocity_change']<gates['phase_velocity_change'] and row['width_relative_change']<gates['width_relative_change'] and
                row['relative_energy_drift']<gates['relative_energy_drift'] and row['zero_axion_exact'] and
                (n!=4096 or row['continuum_field_error']<gates['finest_continuum_field_error']))
            rows.append(row)
    report=dict(kind=__doc__,scope='Initial t=0..4 transient only; not controls over the full packet-scattering window',
        gates=gates,rows=rows,total_solver_seconds=solver_seconds,passed=all(r['passed'] for r in rows),build=build,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_planar_empty_control.py','planar_diagnostics.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
