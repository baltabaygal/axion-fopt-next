"""Moving second-wall control with a predeclared local comparison region."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from planar_boundary_initial import extended_seed,moving_final_core
from planar_diagnostics import wall_geometry
from paper_reference import Parameters

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--suite',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale backend')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    suite=json.loads((args.suite/'report.json').read_text())
    if not all(suite['checks'][key] for key in ('pilot_gates','time_order','space_order','time_size','amplitude')):
        raise RuntimeError('Moving refinement prerequisite failed')
    baseline=json.loads((args.baseline/'report.json').read_text())
    spatial_order=baseline.get('spatial_order',2)
    if suite.get('spatial_order',2)!=spatial_order:raise RuntimeError('Refinement spatial order mismatch')
    if spatial_order not in build.get('supports_line_spatial_orders',[2]):raise RuntimeError('Unsupported spatial order')
    if baseline['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Baseline binary mismatch')
    gates=dict(maximum_total_solver_seconds=15.,initial_core_exact=True,tail_change_energy_relative_max=1e-10,
        final_core_scaled_difference_max=.003,energy_gain_difference_max=2e-4,recoil_difference_max=2e-5,
        continuity_max=2e-3,final_coverage_min=.98)
    final_core=moving_final_core(.3,180.)
    contract=dict(gates=gates,spatial_order=spatial_order,physical_domain=[-256.,512.],original_domain=[0.,512.],
        target_wall=256.,second_wall=-256.,original_second_wall=0.,packet_centre=128.,
        physical_detector_centre=88.,inner=72.,outer=104.,dx=.125,dt=.0125,time_end=180.,speed=.3,initial_common_core=[32.,288.],final_comparison_core=final_core)
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    x,initial,raw=extended_seed(speed=.3);n=len(initial);dx=.125;box=768.;origin=-256.;p=Parameters(dv=0.)
    baseline_root=args.baseline/'packet';old_initial=np.loadtxt(baseline_root/'native/initial.txt')
    old_x=(np.arange(len(old_initial))+.5)*dx;old_mask=(old_x>=32.)&(old_x<=288.);mask=(x>=32.)&(x<=288.)
    initial_exact=bool(np.array_equal(x[mask],old_x[old_mask]) and np.array_equal(initial[mask],old_initial[old_mask]))
    baseline_h=np.genfromtxt(baseline_root/'native/history.csv',delimiter=',',names=True);e0=baseline_h['axion_energy'][0]
    difference=initial[:,[1,3]]-raw
    tail_change=float(np.sum(difference[:,1]**2+(np.roll(difference[:,0],-1)-difference[:,0])**2/dx**2+p.mass_false**2*difference[:,0]**2)*dx/(2*e0))
    if not initial_exact or tail_change>gates['tail_change_energy_relative_max']:raise RuntimeError('Initial overlap/tail gate failed; no solver launched')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    data={};spent=0.
    for label in ('empty','packet'):
        directory=args.output/label;directory.mkdir();s=initial.copy()
        if label=='empty':s[:,1]=0.;s[:,3]=0.
        np.savetxt(directory/'input.txt',s,fmt='%.17g')
        params=dict(line_spatial_order=spatial_order,N=n,lSide=box,dt=.0125,tMax=180.,expansion=False,evolver='PV2',PS_type=2,
            **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,output_every=20,
            line_window_inner=72.,line_window_outer=104.,line_window_centre=88.-origin,line_window_velocity=0.,
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
        remaining=gates['maximum_total_solver_seconds']-spent
        if remaining<=0:raise RuntimeError('Runtime budget exhausted')
        start=time.monotonic();run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],
            env=environment,capture_output=True,text=True,timeout=remaining)
        spent+=time.monotonic()-start;(directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        data[label]=(np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True),np.loadtxt(directory/'native/final.txt'))
        print(label,'completed; solver seconds',spent,flush=True)
    h,final=data['packet'];he,empty=data['empty'];t=h['t'];energy0=h['axion_energy'][0]
    integrate=lambda v:np.r_[0.,np.cumsum(.5*(v[1:]+v[:-1])*np.diff(t))]
    flux=h['window_axion_flux'];q=h['window_exchange']
    residual=h['window_axion_energy']-h['window_axion_energy'][0]+integrate(flux)-integrate(q)
    returned=float(np.trapezoid(np.maximum(-flux,0.),t)/energy0)
    energy_gain=float(h['axion_energy'][-1]/energy0)
    baseline_recoil=baseline['packet_wall_final']['position']-baseline['empty_wall_final']['position']
    old_final=np.loadtxt(baseline_root/'native/final.txt');scales=np.array([1.,.001,1.,.0005])
    final_mask=(x>=final_core[0])&(x<=final_core[1]);old_final_mask=(old_x>=final_core[0])&(old_x<=final_core[1])
    field_difference=float(np.max(abs(final[final_mask]-old_final[old_final_mask])/scales))
    axion_difference=float(np.max(abs(final[final_mask][:,[1,3]]-old_final[old_final_mask][:,[1,3]])/scales[[1,3]]))
    recoil=float(wall_geometry(final,box)['position']-wall_geometry(empty,box)['position'])
    continuity=float(np.max(abs(residual))/energy0);coverage=float(h['window_axion_energy'][-1]/h['axion_energy'][-1])
    checks=dict(initial_core=initial_exact,tail_change=tail_change<gates['tail_change_energy_relative_max'],
        final_core=field_difference<gates['final_core_scaled_difference_max'],energy_gain=abs(energy_gain-baseline['observed_lab_energy_ratio'])<gates['energy_gain_difference_max'],
        recoil=abs(recoil-baseline_recoil)<gates['recoil_difference_max'],continuity=continuity<gates['continuity_max'],
        coverage=coverage>gates['final_coverage_min'],empty_exact=bool(np.all(empty[:,[1,3]]==0.) and np.all(he['axion_energy']==0.)))
    result=dict(kind=__doc__,contract=contract,n=n,initial_core_exact=initial_exact,tail_change_energy_relative=tail_change,
        final_core_scaled_difference=field_difference,final_axion_scaled_difference=axion_difference,
        returned_crossing_fraction=returned,lab_energy_gain=energy_gain,lab_energy_gain_difference=energy_gain-baseline['observed_lab_energy_ratio'],
        recoil=recoil,recoil_difference=recoil-baseline_recoil,continuity_relative=continuity,final_coverage=coverage,
        total_solver_seconds=spent,checks=checks,declared_checks_passed=all(checks.values()),
        scope='Finite-time sensitivity to moving one second wall. Not an all-boundary proof, asymptotic transmission measurement or 3D stability.',build=build,
        baseline_report_sha256=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest(),
        refinement_report_sha256=hashlib.sha256((args.suite/'report.json').read_bytes()).hexdigest(),
        original_waveform_gate_passed=suite['checks']['direct_waveform'],
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('run_moving_boundary_control.py','planar_boundary_initial.py','planar_initial.py','planar_diagnostics.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('build','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
