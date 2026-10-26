"""Bounded time/space controls of the observed planar wall turnaround."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from refine_planar_scattering import coarse_sample

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    base=json.loads((args.baseline/'report.json').read_text())
    if not base['numerical_gates_passed'] or not base['turnaround_observed']:raise RuntimeError('Pressure-pilot prerequisites absent')
    if hashlib.sha256(args.executable.read_bytes()).hexdigest()!=base['build']['executable_sha256']:raise RuntimeError('Baseline executable mismatch')
    bias=float(base['contract'].get('delta_v',0.))
    if base['contract']['n']!=4096 or base['contract']['dt']!=.0125:raise ValueError('Expected central grid/time baseline')
    gates=dict(maximum_wall_seconds=30.,time_ratio_min=3.,space_ratio_min=12.,finest_scaled_field_difference_max=.02,
        finest_trajectory_position_difference_max=.02,finest_turn_time_difference_max=.5)
    cases=[('time-coarse',4096,.025),('time-fine',4096,.00625),('space-coarse',2048,.0125),('space-fine',8192,.0125)]
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,cases=cases,amplitude=.16,
        delta_v=bias,fixed_physics='L=512,t=180,u=.3,potential and continuum packet coefficients unchanged',
        comparison_scales=[1.,.16,1.,.08],classification_gates='Same .02 positive velocity for 5 time units, .5 late position growth'),indent=2)+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    reports={'base':base};states={'base':np.loadtxt(args.baseline/'amplitude-0.16/native/final.txt')}
    hists={'base':np.genfromtxt(args.baseline/'amplitude-0.16/native/history.csv',delimiter=',',names=True)}
    started=time.monotonic()
    for name,n,dt in cases:
        directory=args.output/name;remaining=gates['maximum_wall_seconds']-(time.monotonic()-started)
        if remaining<=0:raise RuntimeError('Pressure refinement runtime budget exhausted')
        command=[sys.executable,str(HERE/'run_planar_pressure_pilot.py'),'--output',str(directory),'--executable',str(args.executable.resolve()),
            '--n',str(n),'--dt',str(dt),'--amplitudes','0','.16','--delta-v',str(bias)]
        run=subprocess.run(command,env=env,capture_output=True,text=True,timeout=remaining)
        (args.output/(name+'.log')).write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError('Case failed: '+name+'; see saved log')
        if not (directory/'report.json').exists():raise RuntimeError('Missing case report: '+name)
        reports[name]=json.loads((directory/'report.json').read_text())
        states[name]=np.loadtxt(directory/'amplitude-0.16/native/final.txt')
        hists[name]=np.genfromtxt(directory/'amplitude-0.16/native/history.csv',delimiter=',',names=True)
        print(name,'completed; wall seconds',time.monotonic()-started,flush=True)
    scales=np.array([1.,.16,1.,.08])
    def difference(a,b):
        if len(a)!=len(b):b=coarse_sample(b,len(a))
        return float(np.max(abs(a-b)/scales))
    td=[difference(states['time-coarse'],states['base']),difference(states['base'],states['time-fine'])]
    sd=[difference(states['space-coarse'],states['base']),difference(states['base'],states['space-fine'])]
    rows={}
    for name,r in reports.items():
        strong=next(v for v in r['rows'] if v['amplitude']==.16)
        rows[name]=dict(n=r['contract']['n'],dt=r['contract']['dt'],final_position=strong['wall_final']['position'],
            final_phase_velocity=strong['wall_final']['phase_velocity'],turn_time=strong['minimum_position_time'],
            minimum_position=strong['minimum_wall_position'],positive_duration=strong['positive_velocity_duration'],
            late_growth=strong['last_position_growth'],observed_turnaround=strong['observed_turnaround'],
            energy_closure=strong['matched_energy_closure'],exchange_balance=strong['exchange_balance'])
    trajectory_difference=float(np.max(abs(hists['base']['wall_x']-hists['space-fine']['wall_x'])))
    turn_difference=abs(rows['base']['turn_time']-rows['space-fine']['turn_time'])
    checks=dict(numerics=all(r['numerical_gates_passed'] for r in reports.values()),
        reversal_persists=all(r['observed_turnaround'] for r in rows.values()),
        time_order=td[0]/td[1]>gates['time_ratio_min'],space_order=sd[0]/sd[1]>gates['space_ratio_min'],
        time_size=td[-1]<gates['finest_scaled_field_difference_max'],space_size=sd[-1]<gates['finest_scaled_field_difference_max'],
        trajectory=trajectory_difference<gates['finest_trajectory_position_difference_max'],timing=turn_difference<gates['finest_turn_time_difference_max'])
    result=dict(kind=__doc__,delta_v=bias,gates=gates,rows=rows,time_differences=td,space_differences=sd,
        time_ratio=td[0]/td[1],space_ratio=sd[0]/sd[1],trajectory_max_difference=trajectory_difference,turn_time_difference=turn_difference,
        checks=checks,declared_checks_passed=all(checks.values()),boundary_control_passed=False,turnaround_validated=False,
        three_dimensional_survival_validated=False,total_solver_seconds=sum(r['total_solver_seconds'] for k,r in reports.items() if k!='base'),
        elapsed_wall_seconds=time.monotonic()-started,baseline_report_sha256=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('refine_planar_pressure.py','run_planar_pressure_pilot.py','refine_planar_scattering.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
