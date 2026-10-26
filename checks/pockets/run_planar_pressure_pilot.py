"""Bounded exploratory field-pressure test; no stable-pocket claim."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from planar_initial import PlanarSeed,prepare
from planar_diagnostics import wall_geometry
from planar_biased_initial import prepare_biased
from paper_reference import Parameters,terms
from planar_stencils import gradient_density

HERE=Path(__file__).resolve().parent


def sustained_interval(mask,t):
    starts=np.flatnonzero(mask&~np.r_[False,mask[:-1]])
    ends=np.flatnonzero(mask&~np.r_[mask[1:],False])
    if len(starts)==0:return 0.,None
    lengths=t[ends]-t[starts];j=int(np.argmax(lengths))
    return float(lengths[j]),[float(t[starts[j]]),float(t[ends[j]])]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--n',type=int,default=4096)
    parser.add_argument('--delta-v',type=float,default=0.)
    parser.add_argument('--dt',type=float,default=.0125)
    parser.add_argument('--amplitudes',type=float,nargs='+',default=[0.,.08,.16])
    args=parser.parse_args()
    if not np.isfinite(args.delta_v) or args.delta_v<0:raise ValueError('Require finite nonnegative bias')
    parameters=Parameters(dv=args.delta_v)
    args.output.mkdir(parents=True,exist_ok=False)
    exe=args.executable.resolve();build=json.loads((exe.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale pressure build')
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Binary mismatch')
    if not build.get('supports_line_wall_geometry') or 4 not in build.get('supports_line_spatial_orders',[]):raise RuntimeError('Required diagnostics/operator unavailable')
    if (len(args.amplitudes)<2 or args.amplitudes[0]!=0. or not np.isfinite(args.amplitudes).all() or
        args.amplitudes!=sorted(set(args.amplitudes)) or not np.isfinite(args.dt) or args.dt<=0 or
        abs(.25/args.dt-round(.25/args.dt))>1e-8):raise ValueError('Require sorted amplitudes starting at zero and integral output cadence')
    gates=dict(maximum_solver_seconds=15.,energy_closure_relative=.002,exchange_balance_relative=.002,
        geometry_endpoint_absolute=1e-12,positive_velocity_threshold=.02,positive_duration_min=5.,
        last_position_interval=20.,last_position_growth_min=.5,moderate_slowing_min=.05)
    contract=dict(gates=gates,n=args.n,box=512.,dt=args.dt,duration=180.,speed=.3,spatial_order=4,
        amplitudes=args.amplitudes,delta_v=args.delta_v,output_interval=.25,
        initial_profile='Boosted degenerate profile; approximate when delta_v>0',
        hypothesis='Stronger axion field pressure can slow and reverse the target wall',
        scope='Planar initial-velocity/bias wall pilot, not long-lived pocket or 3D equilibrium',
        reference_policy='Constant-speed weak-wave Doppler/coverage predictions are not success gates for nonlinear pressure')
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    histories={};rows=[];spent=0.;n=args.n;box=512.
    for amplitude in contract['amplitudes']:
        label='empty' if amplitude==0 else f'amplitude-{amplitude}'
        directory=args.output/label;directory.mkdir()
        seed=PlanarSeed(speed=.3,amplitude=amplitude)
        x,s=prepare_biased(n,seed,parameters) if args.delta_v>0 else prepare(n,seed,parameters)
        np.savetxt(directory/'input.txt',s,fmt='%.17g')
        params=dict(N=n,lSide=box,dt=args.dt,tMax=180.,expansion=False,evolver='PV2',PS_type=2,line_spatial_order=4,
            line_wall_diagnostics=True,**{'lambda':8.},delta_v=args.delta_v,mass_fv=.1,mass_kappa=float(np.log(100.)),axion_fa=10.,output_every=round(.25/args.dt),
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
        remaining=gates['maximum_solver_seconds']-spent
        if remaining<=0:raise RuntimeError('Pressure runtime budget exhausted')
        start=time.monotonic();run=subprocess.run([str(exe),'input='+str((directory/'input.in').resolve())],env=env,
            capture_output=True,text=True,timeout=remaining)
        elapsed=time.monotonic()-start;spent+=elapsed;(directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True);final=np.loadtxt(directory/'native/final.txt')
        before=wall_geometry(s,box);after=wall_geometry(final,box)
        geometry_error=max(abs(h['wall_x'][0]-before['position']),abs(h['wall_x'][-1]-after['position']),
            abs(h['wall_phase_velocity'][-1]-after['phase_velocity']),abs(h['wall_width'][-1]-after['width_10_90']))
        t=h['t'];position=np.unwrap(h['wall_x']*2*np.pi/box)*box/(2*np.pi)
        duration,interval=sustained_interval(h['wall_phase_velocity']>gates['positive_velocity_threshold'],t)
        start_late=int(np.searchsorted(t,t[-1]-gates['last_position_interval']))
        late_growth=float(position[-1]-position[start_late])
        reversal=duration>=gates['positive_duration_min'] and late_growth>gates['last_position_growth_min']
        row=dict(amplitude=amplitude,elapsed_seconds=elapsed,initial_axion_energy=float(h['axion_energy'][0]),
            final_axion_energy=float(h['axion_energy'][-1]),wall_initial=before,wall_final=after,
            minimum_wall_position=float(position.min()),minimum_position_time=float(t[np.argmin(position)]),
            positive_velocity_duration=duration,positive_velocity_interval=interval,last_position_growth=late_growth,
            observed_turnaround=bool(reversal),geometry_endpoint_error=geometry_error,
            minimum_wall_width_cells=float(h['wall_width_cells'].min()))
        if amplitude==0.:
            row['empty_axion_exact']=bool(np.all(final[:,[1,3]]==0.) and np.all(h['axion_energy']==0.))
        else:
            empty=histories[0.];e0=h['axion_energy'][0];q=h['exchange']
            integral=np.r_[0.,np.cumsum(.5*(q[1:]+q[:-1])*np.diff(t))]
            matched=(h['phi_energy']-h['phi_energy'][0])-(empty['phi_energy']-empty['phi_energy'][0])
            row.update(final_axion_ratio=float(h['axion_energy'][-1]/e0),
                matched_energy_closure=float(np.max(abs(matched+h['axion_energy']-e0))/e0),
                exchange_balance=float(np.max(abs(h['axion_energy']-e0-integral))/e0),
                late_mean_velocity_difference_from_empty=float(np.mean(h['wall_phase_velocity'][start_late:]-empty['wall_phase_velocity'][start_late:])))
        density=.5*final[:,3]**2+terms(final[:,0],final[:,1],parameters)[1]+gradient_density(final[:,1],box/n,4)
        row['final_axion_energy_in_tv_region']=float(np.sum(density[final[:,0]>.9])*(box/n))
        row['tv_region_policy']='Endpoint phi>.9; includes penetration and signed gradient split, not cumulative leakage'
        row['mean_contour_acceleration']=(float(h['wall_phase_velocity'][-1])-float(h['wall_phase_velocity'][0]))/float(t[-1]-t[0])
        histories[amplitude]=h;rows.append(row);print(label,'completed; solver seconds',spent,flush=True)
    numerical=dict(empty_exact=rows[0]['empty_axion_exact'],geometry=all(r['geometry_endpoint_error']<gates['geometry_endpoint_absolute'] for r in rows),
        energy=all(r['matched_energy_closure']<gates['energy_closure_relative'] for r in rows[1:]),
        exchange=all(r['exchange_balance']<gates['exchange_balance_relative'] for r in rows[1:]),wall_resolution=all(r['minimum_wall_width_cells']>=4 for r in rows))
    hypothesis=dict(stronger_reverses=rows[-1]['observed_turnaround'],empty_does_not_reverse=not rows[0]['observed_turnaround'])
    moderate=next((r for r in rows if r['amplitude']==.08),None)
    if moderate is not None:hypothesis['moderate_slows']=moderate['late_mean_velocity_difference_from_empty']>gates['moderate_slowing_min']
    result=dict(kind=__doc__,contract=contract,rows=rows,total_solver_seconds=spent,numerical_checks=numerical,
        numerical_gates_passed=all(numerical.values()),hypothesis_checks=hypothesis,
        turnaround_observed=all(hypothesis.values()),turnaround_validated=False,three_dimensional_survival_validated=False,
        remaining=['Separate spatial/time refinement','Amplitude trend beyond two values','Displaced-second-wall control','Distributed trapped fields before spherical pockets'],build=build,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_planar_pressure_pilot.py','planar_initial.py','planar_biased_initial.py','paper_reference.py','planar_stencils.py','planar_diagnostics.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('build','source_sha256')},indent=2))
    if not all(numerical.values()):raise SystemExit(1)

if __name__=='__main__':main()
