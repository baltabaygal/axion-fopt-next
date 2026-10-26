"""Bounded native coupled-field test of moving-wall Doppler and energy gain."""
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
from moving_packet_reference import left_asymptotic_fields

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--n',type=int)
    parser.add_argument('--dt',type=float)
    parser.add_argument('--amplitude',type=float,default=.001)
    parser.add_argument('--spatial-order',type=int,choices=(2,4),default=2)
    parser.add_argument('--require-waveform',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    ref=json.loads(args.reference.read_text())
    if not ref['reference_design_passed'] or not all(ref.get('checks',{}).get(key,False) for key in ('initial_energy','initial_core')):
        raise RuntimeError('Reference design failed initial-packet purity or other checks')
    c=dict(ref['contract'])
    if args.n is not None:c['n']=args.n
    if args.dt is not None:c['dt']=args.dt
    n=c['n'];box=c['box'];dt=c['dt'];duration=c['time_end'];window=c['detector']
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale backend')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    if args.spatial_order not in build.get('supports_line_spatial_orders',[2]):raise RuntimeError('Unsupported compiled spatial order')
    if not np.isfinite(args.amplitude) or args.amplitude<=0:raise ValueError('Positive packet amplitude required')
    gates=dict(direct_waveform_error=.01,maximum_solver_seconds=15.,energy_gain_relative_error=.01,frequency_relative_error=.01,
        axion_continuity_relative=.002,matched_sector_energy_closure=.002,final_coverage_min=.98)
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,physical=c,amplitude=args.amplitude,spatial_order=args.spatial_order,require_waveform=args.require_waveform,reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest()),indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    p=Parameters(dv=0.);seed=PlanarSeed(speed=c['speed'],amplitude=args.amplitude);x,initial=prepare(n,seed,p);data={};spent=0.
    for label in ('empty','packet'):
        directory=args.output/label;directory.mkdir();s=initial.copy()
        if label=='empty':s[:,1]=0.;s[:,3]=0.
        np.savetxt(directory/'input.txt',s,fmt='%.17g')
        params=dict(line_spatial_order=args.spatial_order,N=n,lSide=box,dt=dt,tMax=duration,expansion=False,evolver='PV2',PS_type=2,
            **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,
            output_every=round(.25/dt),line_window_inner=window['inner'],line_window_outer=window['outer'],
            line_window_centre=window['centre'],line_window_velocity=0.,
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
        remaining=gates['maximum_solver_seconds']-spent
        if remaining<=0:raise RuntimeError('Runtime budget exhausted')
        start=time.monotonic();run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],
            env=environment,capture_output=True,text=True,timeout=remaining)
        spent+=time.monotonic()-start;(directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        data[label]=(np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True),np.loadtxt(directory/'native/final.txt'))
        print(label,'completed; solver seconds',spent,flush=True)
    h,s=data['packet'];he,se=data['empty'];e0=h['axion_energy'][0];gain=float(h['axion_energy'][-1]/e0)
    k=2*np.pi*np.fft.fftfreq(n,d=box/n);omega=np.sqrt(p.mass_false**2+k*k)
    ak=np.fft.fft(s[:,1])/n;pk=np.fft.fft(s[:,3])/n
    spectral_energy=.5*(abs(pk)**2+omega**2*abs(ak)**2)
    mean=float(np.sum(omega*spectral_energy)/spectral_energy.sum())
    t=h['t'];integrate=lambda v:np.r_[0.,np.cumsum(.5*(v[1:]+v[:-1])*np.diff(t))]
    residual=h['window_axion_energy']-h['window_axion_energy'][0]+integrate(h['window_axion_flux'])-integrate(h['window_exchange'])
    matched_phi=(h['phi_energy'][-1]-h['phi_energy'][0])-(he['phi_energy'][-1]-he['phi_energy'][0])
    closure=float(abs(matched_phi+h['axion_energy'][-1]-e0)/e0)
    energy_error=abs(gain/ref['reflected_lab_energy_ratio']-1)
    frequency_error=abs(mean/ref['reflected_energy_weighted_frequency']-1)
    coverage=float(h['window_axion_energy'][-1]/h['axion_energy'][-1]);continuity=float(np.max(abs(residual))/e0)
    checks=dict(energy_gain=energy_error<gates['energy_gain_relative_error'],frequency=frequency_error<gates['frequency_relative_error'],
        continuity=continuity<gates['axion_continuity_relative'],sector_closure=closure<gates['matched_sector_energy_closure'],
        coverage=coverage>gates['final_coverage_min'],empty_exact=bool(np.all(se[:,[1,3]]==0.) and np.all(he['axion_energy']==0.)))
    core=(x>=32)&(x<=184)
    reference_a,reference_pa=left_asymptotic_fields(x[core],duration,ref['mode_predictions'],seed)
    field_error=max(float(np.max(abs(s[core,1]-reference_a))/seed.amplitude),
        float(np.max(abs(s[core,3]-reference_pa)))/(seed.omega*seed.amplitude))
    if args.require_waveform:checks['direct_waveform']=field_error<gates['direct_waveform_error']
    report=dict(kind=__doc__,spatial_order=args.spatial_order,scope='Single-grid moving-wall pilot, refinements and boundary/amplitude controls remain pending',
        n=n,box=box,dt=dt,duration=duration,speed=seed.speed,amplitude=seed.amplitude,total_solver_seconds=spent,
        observed_lab_energy_ratio=gain,reference_lab_energy_ratio=ref['reflected_lab_energy_ratio'],energy_gain_relative_error=energy_error,
        observed_energy_weighted_frequency=mean,reference_energy_weighted_frequency=ref['reflected_energy_weighted_frequency'],frequency_relative_error=frequency_error,
        matched_phi_energy_change_over_incident=float(matched_phi/e0),matched_sector_energy_closure=closure,
        final_window_coverage=coverage,axion_continuity_relative=continuity,scaled_reference_field_error=field_error,
        packet_wall_final=wall_geometry(s,box),empty_wall_final=wall_geometry(se,box),
        gates=gates,checks=checks,declared_checks_passed=all(checks.values()),moving_scattering_validated=False,build=build,
        reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('run_moving_packet_pilot.py','moving_packet_reference.py','planar_initial.py','paper_reference.py','planar_diagnostics.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('build','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
