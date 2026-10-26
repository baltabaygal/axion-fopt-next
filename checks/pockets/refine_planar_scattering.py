"""Bounded time/space/amplitude controls for the stationary planar packet.

Boundary variation is a separate missing requirement; this suite cannot
certify isolated scattering or any three-dimensional pocket property.
"""
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
from paper_reference import Parameters

HERE=Path(__file__).resolve().parent


def coarse_sample(fine,n,box=512.):
    """Periodic Fourier shift, then subsample at the coarse cell centres."""
    nf=len(fine)
    if nf%n:raise ValueError('Require integer grid refinement')
    k=2*np.pi*np.fft.fftfreq(nf,d=box/nf)
    delta=.5*(box/n-box/nf)
    shifted=np.fft.ifft(np.fft.fft(fine,axis=0)*np.exp(1j*k[:,None]*delta),axis=0).real
    return shifted[::nf//n]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale backend')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    gates=dict(maximum_total_solver_seconds=30.,time_refinement_ratio_min=3.,space_refinement_ratio_min=3.,
        time_finest_difference_max=.002,space_finest_difference_max=.03,amplitude_scaled_difference_max=.003,
        continuity_relative_max=2e-3,cadence_relative_max=2e-3,final_coverage_min=.98)
    cases=[('time-half',2048,.0125,1.),('time-quarter',2048,.00625,1.),
           ('space-mid',4096,.0125,1.),('space-fine',8192,.0125,1.),('amplitude-half',4096,.0125,.5)]
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,cases=cases,
        fixed_physics='L=512,t=220,degenerate lambda=8,mFV=.1,mTV=1,fa=10,packet omega=.5,width=16,A=.001',
        detector=dict(centre=104.,inner=88.,outer=120.),output_interval=.5,
        state_scales=[1.,.001,1.,.0005],boundary_control_pending=True),indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    spent=0.;all_data={};rows=[];p=Parameters(dv=0.)
    def load(native):return np.genfromtxt(native/'history.csv',delimiter=',',names=True),np.loadtxt(native/'final.txt')
    baseline_report=json.loads((args.baseline/'report.json').read_text())
    if baseline_report['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Baseline binary mismatch')
    all_data['base']=dict(packet=load(args.baseline/'packet/native'),empty=load(args.baseline/'empty-wall/native'))
    for name,n,dt,amplitude in cases:
        data={};row=dict(name=name,n=n,dt=dt,amplitude_factor=amplitude)
        if amplitude!=1.:data['empty']=all_data['space-mid']['empty']
        for label in ('empty','packet'):
            if label in data:continue
            directory=args.output/(name+'-'+label);directory.mkdir()
            _,state=prepare(n,PlanarSeed(amplitude=.001*amplitude),p)
            if label=='empty':state[:,1]=0.;state[:,3]=0.
            np.savetxt(directory/'input.txt',state,fmt='%.17g')
            params=dict(N=n,lSide=512.,dt=dt,tMax=220.,expansion=False,evolver='PV2',PS_type=2,
                **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,
                output_every=round(.5/dt),line_window_inner=88.,line_window_outer=120.,line_window_centre=104.,line_window_velocity=0.,
                initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
            (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
            remaining=gates['maximum_total_solver_seconds']-spent
            if remaining<=0:raise RuntimeError('Runtime budget exhausted; no further solver launched')
            start=time.monotonic();run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],
                env=environment,capture_output=True,text=True,timeout=remaining)
            elapsed=time.monotonic()-start;spent+=elapsed;(directory/'run.log').write_text(run.stdout+run.stderr)
            if run.returncode:raise RuntimeError(run.stderr[-2000:])
            data[label]=load(directory/'native')
            if not np.array_equal(state,np.loadtxt(directory/'native/initial.txt')):raise RuntimeError('Initial transfer mismatch')
        h,s=data['packet'];he,se=data['empty'];t=h['t'];e0=h['axion_energy'][0]
        integrate=lambda v:np.r_[0.,np.cumsum(.5*(v[1:]+v[:-1])*np.diff(t))]
        flux=h['window_axion_flux'];q=h['window_exchange']
        residual=h['window_axion_energy']-h['window_axion_energy'][0]+integrate(flux)-integrate(q)
        ix=np.arange(0,len(t),2)
        if ix[-1]!=len(t)-1:ix=np.r_[ix,len(t)-1]
        row.update(returned_fraction=float(np.trapezoid(np.maximum(-flux,0),t)/e0),
            outward_fraction=float(np.trapezoid(np.maximum(flux,0),t)/e0),
            final_coverage=float(h['window_axion_energy'][-1]/h['axion_energy'][-1]),
            continuity_relative=float(np.max(abs(residual))/e0),
            cadence_relative=float(max(abs(np.trapezoid(v,t)-np.trapezoid(v[ix],t[ix])) for v in (flux,q))/e0),
            wall_recoil=float(wall_geometry(s,512.)['position']-wall_geometry(se,512.)['position']),
            empty_axion_exact=bool(np.all(se[:,[1,3]]==0.) and np.all(he['axion_energy']==0.)))
        rows.append(row);all_data[name]=data
        print(name,'completed; cumulative solver seconds',spent,flush=True)
    scales=np.array([1.,.001,1.,.0005])
    def difference(a,b):
        if len(a)!=len(b):b=coarse_sample(b,len(a))
        return float(np.max(abs(a-b)/scales))
    state=lambda name:all_data[name]['packet'][1]
    time_diffs=[difference(state('base'),state('time-half')),difference(state('time-half'),state('time-quarter'))]
    space_diffs=[difference(state('time-half'),state('space-mid')),difference(state('space-mid'),state('space-fine'))]
    amplitude_diff=float(np.max(abs(state('space-mid')[:,[1,3]]-2*state('amplitude-half')[:,[1,3]])/scales[[1,3]]))
    comparisons=dict(time_differences=time_diffs,time_ratio=time_diffs[0]/time_diffs[1],
        space_differences=space_diffs,space_ratio=space_diffs[0]/space_diffs[1],amplitude_scaled_difference=amplitude_diff)
    checks=dict(time_order=comparisons['time_ratio']>gates['time_refinement_ratio_min'],
        space_order=comparisons['space_ratio']>gates['space_refinement_ratio_min'],
        time_size=time_diffs[-1]<gates['time_finest_difference_max'],space_size=space_diffs[-1]<gates['space_finest_difference_max'],
        linear_amplitude=amplitude_diff<gates['amplitude_scaled_difference_max'],
        continuity=all(r['continuity_relative']<gates['continuity_relative_max'] for r in rows),
        cadence=all(r['cadence_relative']<gates['cadence_relative_max'] for r in rows),
        coverage=all(r['final_coverage']>gates['final_coverage_min'] for r in rows),
        zero_control=all(r['empty_axion_exact'] for r in rows))
    report=dict(kind=__doc__,gates=gates,rows=rows,comparisons=comparisons,checks=checks,
        declared_checks_passed=all(checks.values()),boundary_control_passed=False,physical_scattering_validated=False,
        total_solver_seconds=spent,build=build,baseline_sha256=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('refine_planar_scattering.py','planar_initial.py','planar_diagnostics.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(rows=rows,comparisons=comparisons,checks=checks,total_solver_seconds=spent),indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
