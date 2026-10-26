"""Tiny compiled formation evolution vs independent semi-discrete dynamics.

Fixed spatial stencil; this establishes temporal/diagnostic correctness, not
spatial convergence, physical trapping, stability or long-lived survival.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from model import config,prepare_run,read_state,radius,laplacian
from paper_reference import Parameters,force,terms
from monitor_reference import monitor,continuity_summary
from tracking_monitor_reference import tracking_monitor

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Backend build is stale: '+name)
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:
        raise RuntimeError('Executable differs from build manifest')
    p=Parameters();n=8;box=8.;dx=box/n;duration=.2
    state=np.zeros((n,n,n,4));r=radius(n,box)
    state[...,0]=.5*(1+np.tanh((r-1.7)/.6));state[...,1]=2.*np.exp(-r*r/3.)
    state[...,2]=-.1*(np.roll(state[...,0],-1,0)-np.roll(state[...,0],1,0))+.02*np.exp(-r*r/3.)
    state[...,3]=.15*np.exp(-r*r/3.)
    scales=np.array([1.,2.,1.,.2])
    gates=dict(finest_scaled_trajectory=2e-5,min_refinement_ratio=3.5,
        snapshot_monitor_relative=1e-12,finest_integrated_continuity=1e-5,
        finest_energy_drift=1e-5,periodic_translation_absolute=1e-12,
        output_cadence_difference=1e-5)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def rhs(t,y):
        s=y.reshape(state.shape);fp,fa=force(s[...,0],s[...,1],p)
        return np.stack((s[...,2],s[...,3],laplacian(s[...,0],dx)-fp,
                         laplacian(s[...,1],dx)-fa),axis=-1).ravel()
    ref=solve_ivp(rhs,(0,duration),state.ravel(),method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    if not ref.success:raise RuntimeError(ref.message)
    np.save(args.output/'reference-final.npy',ref.y[:,-1].reshape(state.shape))
    results=[];histories=[];finals=[]
    def execute(s,dt,directory):
        cfg=config(n,box,dt,duration,lam=p.lam,dv=p.dv,mass_tv=1.,every=1,frames=1,
                   monitor_radius=2.7,tracking_inner=2.4,tracking_outer=3.5)
        cfg.update(potential_model='formation',mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa)
        prepare_run(s,cfg,directory)
        run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],
                           env=environment,capture_output=True,text=True,timeout=30)
        (directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        return np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True)
    for dt in (.01,.005,.0025):
        directory=args.output/f'dt-{dt}';h=execute(state,dt,directory);histories.append(h)
        max_error=max_monitor_error=0.;reference_fields=[]
        for j,t in enumerate(h['u']):
            _,s=read_state(directory/f'native/frame_{j:07d}.bin')
            target=ref.sol(t).reshape(state.shape)
            max_error=max(max_error,float(np.max(np.abs(s-target)/scales)))
            phi,a,pp,pa=(s[...,k] for k in range(4));u,va,_,_=terms(phi,a,p)
            energy=.5*(pp*pp+pa*pa)+u+va
            for axis in range(3):
                energy+=((np.roll(phi,-1,axis)-phi)**2+(np.roll(a,-1,axis)-a)**2)/(2*dx*dx)
            max_monitor_error=max(max_monitor_error,float(abs(energy.sum()*dx**3-h['total'][j])/h['total'][0]))
            fixed=monitor(s,box,2.7,p.lam,p.dv,1.,formation=p)
            moving=tracking_monitor(s,box,2.4,3.5,p.lam,p.dv,1.,formation=p)
            for reference in (fixed,moving):
                for key,value in reference.items():
                    if key in h.dtype.names:
                        max_monitor_error=max(max_monitor_error,float(abs(value-h[key][j])/h['total'][0]))
            reference_fields.append(target)
        np.save(directory/'reference.npy',np.stack(reference_fields))
        finals.append(s)
        summaries={prefix:continuity_summary(h,prefix) for prefix in ('monitor','tracking')}
        results.append(dict(dt=dt,scaled_trajectory=max_error,snapshot_monitor_relative=max_monitor_error,
            energy_drift=float(np.max(np.abs(h['total']-h['total'][0]))/h['total'][0]),
            continuity=summaries,max_integrated_continuity=max(max(v['max_continuity_residual'].values()) for v in summaries.values()),
            max_output_cadence_difference=max(max(v['coarser_cadence_difference'].values()) for v in summaries.values()),
            moving_speed=float(np.max(np.sqrt(sum(h['tracking_v'+axis]**2 for axis in 'xyz')))),
            exchange_amplitude=float(np.max(np.abs(h['tracking_exchange'])))))
    shift=(3,-2,1);translated=execute(np.roll(state,shift,axis=(0,1,2)),.0025,args.output/'periodic-shift')
    _,translated_final=read_state(args.output/'periodic-shift/native/final.bin')
    field_shift_error=float(np.max(np.abs(translated_final-np.roll(finals[-1],shift,axis=(0,1,2)))))
    invariant_columns=['total']+['tracking_'+k for k in ('energy','phi_energy','axion_energy','phi_flux','axion_flux','exchange',
        'lab_phi_flux','lab_axion_flux','phi_advection','axion_advection','vx','vy','vz','phase_volume','coherence_x','coherence_y','coherence_z','max_fv_radius')]
    diagnostic_shift_error=max(float(np.max(np.abs(translated[key]-histories[-1][key]))) for key in invariant_columns)
    fine=results[-1]
    checks=dict(trajectory=fine['scaled_trajectory']<gates['finest_scaled_trajectory'],
        time_refinement=all(results[j]['scaled_trajectory']/results[j+1]['scaled_trajectory']>gates['min_refinement_ratio'] for j in (0,1)),
        snapshot_monitor_agreement=all(v['snapshot_monitor_relative']<gates['snapshot_monitor_relative'] for v in results),
        integrated_continuity=fine['max_integrated_continuity']<gates['finest_integrated_continuity'],
        energy_drift=fine['energy_drift']<gates['finest_energy_drift'],
        independent_cadence_check=fine['max_output_cadence_difference']<gates['output_cadence_difference'],
        periodic_translation=max(field_shift_error,diagnostic_shift_error)<gates['periodic_translation_absolute'],
        nontrivial_exchange_and_motion=all(v['exchange_amplitude']>1e-5 and v['moving_speed']>.01 for v in results))
    names=['validate_formation_inhomogeneous.py','paper_reference.py','model.py','monitor_reference.py','tracking_monitor_reference.py','pocket_geometry.py']
    report=dict(kind=__doc__,grid=n,box=box,duration=duration,time='Minkowski cosmic time',state_scales=scales.tolist(),
        reference=dict(method='DOP853',rtol=2e-12,atol=2e-14,nfev=ref.nfev,spatial_operator='same periodic nearest-neighbour stencil; no spatial convergence claim'),
        gates=gates,results=results,translation=dict(shift=shift,field_error=field_shift_error,diagnostic_error=diagnostic_shift_error),
        checks=checks,passed=all(checks.values()),build=build,source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],results=[{k:v for k,v in r.items() if k!='continuity'} for r in results],translation=report['translation'],checks=checks),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
