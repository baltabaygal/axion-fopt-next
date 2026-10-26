"""Displace the second wall, preserving local stationary scattering data."""
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
from paper_reference import Parameters

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--suite',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale backend')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    suite=json.loads((args.suite/'report.json').read_text())
    if not suite['declared_checks_passed']:raise RuntimeError('Baseline refinement failed')
    if suite['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Baseline binary mismatch')
    gates=dict(maximum_total_solver_seconds=15.,initial_core_exact=True,tail_change_energy_relative_max=1e-10,
        final_core_scaled_difference_max=.003,returned_difference_max=2e-4,recoil_difference_max=2e-5,
        continuity_max=2e-3,final_coverage_min=.98)
    contract=dict(gates=gates,physical_domain=[-256.,512.],original_domain=[0.,512.],
        target_wall=256.,second_wall=-256.,original_second_wall=0.,packet_centre=128.,
        physical_detector_centre=104.,inner=88.,outer=120.,dx=.125,dt=.0125,time_end=220.,common_core=[32.,288.])
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    x,initial,raw=extended_seed();n=len(initial);dx=.125;box=768.;origin=-256.;p=Parameters(dv=0.)
    baseline_root=args.suite/'space-mid-packet';old_initial=np.loadtxt(baseline_root/'native/initial.txt')
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
        params=dict(N=n,lSide=box,dt=.0125,tMax=220.,expansion=False,evolver='PV2',PS_type=2,
            **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,output_every=40,
            line_window_inner=88.,line_window_outer=120.,line_window_centre=104.-origin,line_window_velocity=0.,
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
    reference=next(r for r in suite['rows'] if r['name']=='space-mid')
    old_final=np.loadtxt(baseline_root/'native/final.txt');scales=np.array([1.,.001,1.,.0005])
    field_difference=float(np.max(abs(final[mask]-old_final[old_mask])/scales))
    axion_difference=float(np.max(abs(final[mask][:,[1,3]]-old_final[old_mask][:,[1,3]])/scales[[1,3]]))
    recoil=float(wall_geometry(final,box)['position']-wall_geometry(empty,box)['position'])
    continuity=float(np.max(abs(residual))/energy0);coverage=float(h['window_axion_energy'][-1]/h['axion_energy'][-1])
    checks=dict(initial_core=initial_exact,tail_change=tail_change<gates['tail_change_energy_relative_max'],
        final_core=field_difference<gates['final_core_scaled_difference_max'],returned=abs(returned-reference['returned_fraction'])<gates['returned_difference_max'],
        recoil=abs(recoil-reference['wall_recoil'])<gates['recoil_difference_max'],continuity=continuity<gates['continuity_max'],
        coverage=coverage>gates['final_coverage_min'],empty_exact=bool(np.all(empty[:,[1,3]]==0.) and np.all(he['axion_energy']==0.)))
    result=dict(kind=__doc__,contract=contract,n=n,initial_core_exact=initial_exact,tail_change_energy_relative=tail_change,
        final_core_scaled_difference=field_difference,final_axion_scaled_difference=axion_difference,
        returned_fraction=returned,returned_difference=returned-reference['returned_fraction'],
        recoil=recoil,recoil_difference=recoil-reference['wall_recoil'],continuity_relative=continuity,final_coverage=coverage,
        total_solver_seconds=spent,checks=checks,declared_checks_passed=all(checks.values()),
        scope='Finite-time sensitivity to moving one second wall. Not an all-boundary proof, asymptotic transmission measurement or 3D stability.',build=build,
        baseline_report_sha256=hashlib.sha256((args.suite/'report.json').read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('run_planar_boundary_control.py','planar_boundary_initial.py','planar_initial.py','planar_diagnostics.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('build','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
