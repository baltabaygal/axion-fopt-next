"""Bounded moving-wall time, space and amplitude checks; boundary still pending."""
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
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--spatial-order',type=int,choices=(2,4),default=2)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    baseline=json.loads((args.baseline/'report.json').read_text())
    if baseline.get('spatial_order',2)!=args.spatial_order:raise RuntimeError('Baseline spatial order mismatch')
    if not baseline['declared_checks_passed']:raise RuntimeError('Pilot prerequisite failed')
    gates=dict(maximum_total_wall_seconds=35.,time_ratio_min=3.,space_ratio_min=(12. if args.spatial_order==4 else 3.),
        time_finest_difference_max=.002,amplitude_scaled_difference_max=.003,direct_waveform_error_max=.01)
    cases=[('time-coarse',4096,.025,.001),('time-fine',4096,.00625,.001),
           ('space-coarse',2048,.0125,.001),('space-fine',8192,.0125,.001),('amplitude-half',4096,.0125,.0005)]
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,cases=cases,spatial_order=args.spatial_order,physical='L=512,t=180,u=.3,unchanged potential and continuum packet',
        fixed_scales=[1.,.001,1.,.0005],boundary_pending=True),indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    reports={'base':baseline};states={'base':np.loadtxt(args.baseline/'packet/native/final.txt')};started=time.monotonic()
    for name,n,dt,amplitude in cases:
        directory=args.output/name
        remaining=gates['maximum_total_wall_seconds']-(time.monotonic()-started)
        if remaining<=0:raise RuntimeError('Suite budget exhausted')
        command=[sys.executable,str(HERE/'run_moving_packet_pilot.py'),'--output',str(directory),
            '--executable',str(args.executable.resolve()),'--reference',str(args.reference.resolve()),
            '--n',str(n),'--dt',str(dt),'--amplitude',str(amplitude),'--spatial-order',str(args.spatial_order)]
        run=subprocess.run(command,env=environment,capture_output=True,text=True,timeout=remaining)
        (args.output/(name+'.log')).write_text(run.stdout+run.stderr)
        if not (directory/'report.json').exists():raise RuntimeError('Case terminated without report: '+name)
        report=json.loads((directory/'report.json').read_text());reports[name]=report
        states[name]=np.loadtxt(directory/'packet/native/final.txt')
        print(name,'completed; elapsed seconds',time.monotonic()-started,flush=True)
    scale=np.array([1.,.001,1.,.0005])
    def diff(a,b):
        if len(a)!=len(b):b=coarse_sample(b,len(a))
        return float(np.max(abs(a-b)/scale))
    time_diffs=[diff(states['time-coarse'],states['base']),diff(states['base'],states['time-fine'])]
    space_diffs=[diff(states['space-coarse'],states['base']),diff(states['base'],states['space-fine'])]
    amp=float(np.max(abs(states['base'][:,[1,3]]-2*states['amplitude-half'][:,[1,3]])/scale[[1,3]]))
    comparisons=dict(time_differences=time_diffs,time_ratio=time_diffs[0]/time_diffs[1],
        space_differences=space_diffs,space_ratio=space_diffs[0]/space_diffs[1],amplitude_scaled_difference=amp,
        direct_waveform_errors={name:reports[name]['scaled_reference_field_error'] for name in ('space-coarse','base','space-fine')})
    checks=dict(pilot_gates=all(r['declared_checks_passed'] for r in reports.values()),
        time_order=comparisons['time_ratio']>gates['time_ratio_min'],space_order=comparisons['space_ratio']>gates['space_ratio_min'],
        time_size=time_diffs[-1]<gates['time_finest_difference_max'],amplitude=amp<gates['amplitude_scaled_difference_max'],
        direct_waveform=reports['space-fine']['scaled_reference_field_error']<gates['direct_waveform_error_max'])
    summary={name:{k:r[k] for k in ('n','dt','observed_lab_energy_ratio','observed_energy_weighted_frequency',
        'matched_sector_energy_closure','final_window_coverage','scaled_reference_field_error','declared_checks_passed')} for name,r in reports.items()}
    result=dict(kind=__doc__,spatial_order=args.spatial_order,gates=gates,cases=summary,comparisons=comparisons,checks=checks,
        declared_checks_passed=all(checks.values()),moving_scattering_validated=False,boundary_control_passed=False,
        elapsed_wall_seconds=time.monotonic()-started,total_solver_seconds=sum(r['total_solver_seconds'] for name,r in reports.items() if name!='base'),
        baseline_report_sha256=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('refine_moving_scattering.py','run_moving_packet_pilot.py','moving_packet_reference.py','refine_planar_scattering.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256',)},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
