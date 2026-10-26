"""Small one-thread compiled framework validation of the formation model."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from model import config,prepare_run,read_state
from paper_reference import Parameters,force,terms,hessian

HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();p=Parameters();initial=np.zeros((4,4,4,4));initial[...,:]=[0.,2.,0.,.15]
    gates=dict(finest_scaled_trajectory=2e-5,min_refinement_ratio=3.5,spatial_spread=1e-13,
               finest_energy_drift=1e-5,snapshot_energy_relative=1e-12,frequency_bound_absolute=1e-12)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    times=np.linspace(0,4,21)
    def rhs(t,y):
        fp,fa=force(y[0],y[1],p);return [y[2],y[3],-float(fp),-float(fa)]
    ref=solve_ivp(rhs,(0,4),initial[0,0,0],t_eval=times,method='DOP853',rtol=2e-13,atol=2e-15)
    if not ref.success:raise RuntimeError(ref.message)
    np.savetxt(args.output/'reference.csv',np.column_stack((times,ref.y.T)),delimiter=',',header='t,phi,a,phi_dot,a_dot',comments='')
    def invoke(cfg,state,directory,expect_success=True):
        prepare_run(state,cfg,directory)
        result=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],env=environment,capture_output=True,text=True,timeout=30)
        (directory/'run.log').write_text(result.stdout+result.stderr)
        if expect_success and result.returncode:raise RuntimeError(result.stderr[-2000:])
        return result
    results=[];scales=np.array([1.,2.,1.,.2])
    base=config(4,8.,.02,4.,lam=p.lam,dv=p.dv,mass_tv=1.,every=10,frames=10,monitor_radius=2.7)
    base.update(potential_model='formation',mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa)
    for dt in (.02,.01,.005):
        cfg=dict(base,dt=dt,output_every=round(.2/dt),snapshot_every=round(.2/dt))
        directory=args.output/f'dt-{dt}';invoke(cfg,initial,directory)
        native=directory/'native';history=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
        error=spread=energy_error=frequency_error=0.
        np.testing.assert_allclose(history['u'],times,rtol=0.,atol=1e-14)
        for j,t in enumerate(times):
            meta,state=read_state(native/f'frame_{round(t/dt):07d}.bin')
            if abs(meta['u']-t)>1e-14:raise RuntimeError('Wrong snapshot time')
            error=max(error,float(np.max(np.abs(state-ref.y[:,j])/scales)))
            spread=max(spread,float(np.max(np.ptp(state.reshape(-1,4),axis=0))))
            phi,a,pp,pa=state[0,0,0];u,va,_,_=terms(phi,a,p)
            energy=(.5*(pp*pp+pa*pa)+u+va)*8.**3
            energy_error=max(energy_error,float(abs(energy-history['total'][j])/history['total'][0]))
            xx,xy,aa=hessian(phi,a,p)
            largest=.5*(xx+aa+np.hypot(xx-aa,2*xy))
            expected=dt*np.sqrt(12./2.**2+max(0.,largest))
            frequency_error=max(frequency_error,float(abs(expected-history['dt_frequency_bound'][j])))
        _,transferred=read_state(native/'initial.bin')
        results.append(dict(dt=dt,scaled_trajectory=error,spatial_spread=spread,
            snapshot_energy_relative=energy_error,frequency_bound_absolute=frequency_error,
            energy_drift=float(np.max(np.abs(history['total']-history['total'][0]))/history['total'][0]),
            initial_transfer_exact=bool(np.array_equal(initial,transferred))))
    checks=dict(finest_scaled_trajectory=results[-1]['scaled_trajectory']<gates['finest_scaled_trajectory'],
        timestep_refinement=all(results[i]['scaled_trajectory']/results[i+1]['scaled_trajectory']>gates['min_refinement_ratio'] for i in (0,1)),
        homogeneous_spread=all(v['spatial_spread']<gates['spatial_spread'] for v in results),
        snapshot_energy=all(v['snapshot_energy_relative']<gates['snapshot_energy_relative'] for v in results),
        mixed_hessian_frequency=all(v['frequency_bound_absolute']<gates['frequency_bound_absolute'] for v in results),
        finest_energy_drift=results[-1]['energy_drift']<gates['finest_energy_drift'],
        initial_transfer=all(v['initial_transfer_exact'] for v in results))
    rejections=[]
    for name,changes in [('wrong-tv',dict(mass_tv=2.)),('zero-fa',dict(axion_fa=0.)),
                         ('unknown-model',dict(potential_model='bad')),('polynomial-power',dict(mass_power=4)),
                         ('polynomial-stiffening',dict(vacuum_stiffening=12.)),('missing-fv',dict(mass_fv=None))]:
        cfg=dict(base,tMax=.02);cfg.update(changes);cfg={k:v for k,v in cfg.items() if v is not None}
        directory=args.output/('reject-'+name);run=invoke(cfg,initial,directory,False)
        rejections.append(dict(case=name,rejected=run.returncode!=0,no_native_output=not (directory/'native').exists()))
    checks['parameter_rejections']=all(v['rejected'] and v['no_native_output'] for v in rejections)
    source_names=['validate_formation_backend.py','paper_reference.py','model.py']
    report=dict(kind=__doc__,scope='homogeneous only; not pocket survival or spatial convergence',
        grid=4,box=8.,time='Minkowski cosmic time 0..4',initial=[0.,2.,0.,.15],state_scales=scales.tolist(),
        gates=gates,results=results,rejections=rejections,checks=checks,passed=bool(all(checks.values())),
        provenance=dict(executable=str(executable),executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
            source_sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in source_names},
            environment={k:environment[k] for k in ('OMP_NUM_THREADS','OMP_DYNAMIC','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}))
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],results=results,checks=checks),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
