"""Matched field-state and end-to-end timings for optional CPU threading.

Contract: THREADED_BACKEND.md. This is a numerical/backend check, not
evidence for physical pocket survival. The original Serial run is untouched.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from model import EXECUTABLE, config, eigen_state, history, radial_modes, read_state, run


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--serial',type=Path,default=EXECUTABLE)
    p.add_argument('--threaded',type=Path,default=Path('build/pockets-omp/pocket_evolve'))
    p.add_argument('--probe',type=Path,default=Path('build/pockets-omp/pocket_backend_probe'))
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    modes=radial_modes(8.,1.,2.,.2,1.2,14.,h=.01,mass_power=4)
    initial=eigen_state(64,32.,8.,1.,modes,.8,0,14.)
    cfg=config(64,32.,.0125,2.,lam=8.,dv=.01,mass_tv=2.,mass_power=4,
               every=20,monitor_radius=12.)
    contract=dict(config=cfg,state='Bound-mode Gaussian field, R0=8,w=1,seed=0,amplitude=.8',
                   executable_hashes=dict(serial=sha(args.serial),threaded=sha(args.threaded),probe=sha(args.probe)),
                   repeats=3,threads=[1,2,4],state_tolerance=1e-10,history_tolerance=1e-10,
                   source_sha256=sha(__file__),timing='End-to-end including identical input/output work; concurrent Serial physics run may affect times.')
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    before=time.perf_counter();run(initial,cfg,args.output/'serial',args.serial.resolve())
    serial_seconds=time.perf_counter()-before
    _,serial_state=read_state(args.output/'serial/native/final.bin');serial_history=history(args.output/'serial')
    state_scale=max(1.,float(np.max(np.abs(serial_state))))
    energy_scale=abs(float(serial_history['total'][0]));results={}
    for threads in (1,2,4):
        os.environ['OMP_NUM_THREADS']=str(threads);os.environ['OMP_PROC_BIND']='false'
        probe=subprocess.run([str(args.probe.resolve())],capture_output=True,text=True,check=True)
        (args.output/f'probe_{threads}.log').write_text(probe.stdout+probe.stderr)
        lines=[line.removeprefix('POCKET_BACKEND ') for line in probe.stdout.splitlines() if line.startswith('POCKET_BACKEND ')]
        actual=json.loads(lines[-1])
        if actual['name']!='OpenMP' or actual['concurrency']!=threads:
            raise RuntimeError(f'Requested thread count not active: {actual}')
        samples=[];errors=[]
        for repeat in range(3):
            destination=args.output/f'omp_{threads}_{repeat}'
            before=time.perf_counter();run(initial,cfg,destination,args.threaded.resolve())
            elapsed=time.perf_counter()-before
            _,state=read_state(destination/'native/final.bin');h=history(destination)
            state_error=float(np.max(np.abs(state-serial_state))/state_scale)
            history_error=max(float(np.max(np.abs(h[key]-serial_history[key])))/energy_scale
                              for key in h.dtype.names if key!='u')
            samples.append(elapsed);errors.append(dict(state=state_error,history=history_error))
        results[str(threads)]=dict(actual=actual,seconds=samples,median_seconds=float(np.median(samples)),
                                   agreement=errors,passed=all(e['state']<1e-10 and e['history']<1e-10 for e in errors))
    report=dict(results=results,serial_seconds=serial_seconds,
                agreement_passed=all(row['passed'] for row in results.values()),
                preferred_threads=min(results,key=lambda key:results[key]['median_seconds']),
                interpretation='Matched modest-state timing only. Do not assume the same speedup on larger grids.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
