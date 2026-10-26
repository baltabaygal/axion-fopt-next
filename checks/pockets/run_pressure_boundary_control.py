"""Displaced-second-wall control for the finite-amplitude planar reversal."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from planar_boundary_initial import extended_seed
from planar_diagnostics import wall_geometry
from run_planar_pressure_pilot import sustained_interval

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--refinement',type=Path,required=True);parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    base=json.loads((args.baseline/'report.json').read_text());refined=json.loads((args.refinement/'report.json').read_text())
    if not base['numerical_gates_passed'] or not refined['declared_checks_passed']:raise RuntimeError('Pressure prerequisites failed')
    bias=float(base['contract'].get('delta_v',0.))
    if float(refined.get('delta_v',0.))!=bias:raise RuntimeError('Refinement bias mismatch')
    if refined['baseline_report_sha256']!=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest():raise RuntimeError('Refinement baseline mismatch')
    core_low=80. if bias>0 else 64.
    exe=args.executable.resolve();build=json.loads((exe.parent/'build.json').read_text())
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=base['build']['executable_sha256']:raise RuntimeError('Baseline binary mismatch')
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale pressure build')
    gates=dict(maximum_solver_seconds=15.,initial_core_exact=True,local_state_difference_max=.003,
        wall_trajectory_difference_max=.002,turn_time_difference_max=.5,energy_ratio_difference_max=2e-4,
        matched_energy_closure_max=.002,exchange_balance_max=.002)
    contract=dict(gates=gates,initial_common_core=[32.,288.],final_comparison_core=[core_low,288.],delta_v=bias,
        region_policy='Exclude moved partner: biased core starts at 80, beyond baseline empty partner endpoint ~62.2 plus >16 widths; degenerate core starts at 64',
        physical_domain=[-256.,512.],n=6144,dx=.125,dt=.0125,duration=180.,amplitude=.16,
        reversal_criteria=base['contract']['gates'])
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    x,initial,raw=extended_seed(speed=.3,amplitude=.16);n=len(initial);box=768.
    old_initial=np.loadtxt(args.baseline/'amplitude-0.16/native/initial.txt');xo=(np.arange(len(old_initial))+.5)*.125
    im=(x>=32)&(x<=288);oi=(xo>=32)&(xo<=288)
    exact=bool(np.array_equal(initial[im],old_initial[oi]))
    if not exact:raise RuntimeError('Initial common core differs; no solver launched')
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    data={};spent=0.
    for label in ('empty','packet'):
        directory=args.output/label;directory.mkdir();s=initial.copy()
        if label=='empty':s[:,1]=0.;s[:,3]=0.
        np.savetxt(directory/'input.txt',s,fmt='%.17g')
        params=dict(N=n,lSide=box,dt=.0125,tMax=180.,expansion=False,evolver='PV2',PS_type=2,line_spatial_order=4,line_wall_diagnostics=True,
            **{'lambda':8.},delta_v=bias,mass_fv=.1,mass_kappa=float(np.log(100.)),axion_fa=10.,output_every=20,
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
        remaining=gates['maximum_solver_seconds']-spent
        if remaining<=0:raise RuntimeError('Pressure boundary runtime cap exhausted')
        start=time.monotonic();run=subprocess.run([str(exe),'input='+str((directory/'input.in').resolve())],env=env,capture_output=True,text=True,timeout=remaining)
        spent+=time.monotonic()-start;(directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        data[label]=(np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True),np.loadtxt(directory/'native/final.txt'))
        print(label,'completed; solver seconds',spent,flush=True)
    h,s=data['packet'];he,se=data['empty'];oldh=np.genfromtxt(args.baseline/'amplitude-0.16/native/history.csv',delimiter=',',names=True)
    old=np.loadtxt(args.baseline/'amplitude-0.16/native/final.txt');oldrow=next(r for r in base['rows'] if r['amplitude']==.16)
    mask=(x>=core_low)&(x<=288);omask=(xo>=core_low)&(xo<=288);scales=np.array([1.,.16,1.,.08])
    local=float(np.max(abs(s[mask]-old[omask])/scales));trajectory=float(np.max(abs((h['wall_x']-256.)-oldh['wall_x'])))
    t=h['t'];position=h['wall_x']-256.;turn=float(t[np.argmin(position)])
    duration,interval=sustained_interval(h['wall_phase_velocity']>.02,t)
    growth=float(position[-1]-position[np.searchsorted(t,t[-1]-20.)]);reversal=duration>=5. and growth>.5
    e0=h['axion_energy'][0];matched=(h['phi_energy']-h['phi_energy'][0])-(he['phi_energy']-he['phi_energy'][0])
    closure=float(np.max(abs(matched+h['axion_energy']-e0))/e0);q=h['exchange']
    integral=np.r_[0.,np.cumsum(.5*(q[1:]+q[:-1])*np.diff(t))]
    exchange=float(np.max(abs(h['axion_energy']-e0-integral))/e0);ratio=float(h['axion_energy'][-1]/e0)
    checks=dict(initial=exact,local_fields=local<gates['local_state_difference_max'],trajectory=trajectory<gates['wall_trajectory_difference_max'],
        timing=abs(turn-oldrow['minimum_position_time'])<gates['turn_time_difference_max'],reversal=bool(reversal),
        energy_ratio=abs(ratio-oldrow['final_axion_ratio'])<gates['energy_ratio_difference_max'],energy=closure<gates['matched_energy_closure_max'],
        exchange=exchange<gates['exchange_balance_max'],empty_exact=bool(np.all(se[:,[1,3]]==0.) and np.all(he['axion_energy']==0.)))
    result=dict(kind=__doc__,contract=contract,checks=checks,declared_checks_passed=all(checks.values()),
        local_scaled_field_difference=local,wall_trajectory_max_difference=trajectory,turn_time=turn,
        final_phase_velocity=float(h['wall_phase_velocity'][-1]),positive_duration=duration,positive_interval=interval,late_growth=growth,
        initial_axion_energy=float(e0),final_axion_ratio=ratio,energy_ratio_difference=ratio-oldrow['final_axion_ratio'],
        matched_energy_closure=closure,exchange_balance=exchange,total_solver_seconds=spent,
        scope='One displaced-wall finite-time sensitivity test; not all-boundary or 3D lifetime validation',
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('run_pressure_boundary_control.py','planar_boundary_initial.py','run_planar_pressure_pilot.py')},
        baseline_sha256=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest(),refinement_sha256=hashlib.sha256((args.refinement/'report.json').read_bytes()).hexdigest(),build=build)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('build','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
