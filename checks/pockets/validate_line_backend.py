"""Small correctness gates for actual one-dimensional CosmoLattice fields."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from paper_reference import Parameters,force,terms,hessian
from planar_stencils import laplacian as spatial_laplacian,gradient_density,symbol

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--spatial-order',type=int,choices=(2,4),default=2)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale line build: '+name)
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong executable hash')
    if args.spatial_order not in build.get('supports_line_spatial_orders',[2]):raise RuntimeError('Unsupported compiled spatial order')
    gates=dict(frequency_bound_absolute=1e-12,finest_state_error=2e-5,min_refinement_ratio=3.5,finest_energy_drift=1e-5,
               finest_exchange_residual=1e-5,snapshot_energy_relative=1e-12,translation_absolute=1e-12)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    p=Parameters();n=32;box=8.;dx=box/n;duration=.2;x=(np.arange(n)+.5)*dx;k=2*np.pi/box
    state=np.column_stack((.5+.4*np.sin(k*x),.3*np.cos(k*x)+.2*np.sin(2*k*x),.05*np.cos(k*x),-.1*np.sin(k*x)))
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def rhs(t,y):
        s=y.reshape(n,4);fp,fa=force(s[:,0],s[:,1],p)
        lap=(np.roll(s[:,:2],-1,axis=0)+np.roll(s[:,:2],1,axis=0)-2*s[:,:2])/dx**2
        if args.spatial_order==4:lap=np.column_stack((spatial_laplacian(s[:,0],dx,4),spatial_laplacian(s[:,1],dx,4)))
        return np.column_stack((s[:,2:],lap-np.column_stack((fp,fa)))).ravel()
    reference=solve_ivp(rhs,(0,duration),state.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
    if not reference.success:raise RuntimeError(reference.message)
    ref=reference.y[:,-1].reshape(n,4);np.save(args.output/'reference-final.npy',ref)
    def execute(initial,dt,directory):
        directory.mkdir();np.savetxt(directory/'input.txt',initial,fmt='%.17g')
        params=dict(line_spatial_order=args.spatial_order,N=n,lSide=box,dt=dt,tMax=duration,expansion=False,evolver='PV2',PS_type=2,
            **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,output_every=1,
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{key} = {str(v).lower() if isinstance(v,bool) else v}\n' for key,v in params.items()))
        run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],env=environment,capture_output=True,text=True,timeout=30)
        (directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        return np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True),np.loadtxt(directory/'native/final.txt')
    results=[];histories=[];finals=[]
    for dt in (.01,.005,.0025):
        directory=args.output/f'dt-{dt}';h,s=execute(state,dt,directory);histories.append(h);finals.append(s)
        u,va,_,_=terms(s[:,0],s[:,1],p)
        gradient=(np.roll(s[:,:2],-1,axis=0)-s[:,:2])/dx
        energy=float(np.sum(.5*(s[:,2]**2+s[:,3]**2)+u+va+.5*np.sum(gradient**2,axis=1))*dx)
        if args.spatial_order==4:
            energy=float(np.sum(.5*(s[:,2]**2+s[:,3]**2)+u+va+gradient_density(s[:,0],dx,4)+gradient_density(s[:,1],dx,4))*dx)
        xx,xy,aa=hessian(s[:,0],s[:,1],p)
        maximum=max(0.,float(np.max(.5*(xx+aa+np.hypot(xx-aa,2*xy)))))
        expected_frequency=dt*np.sqrt(symbol(np.pi/dx,dx,args.spatial_order)+maximum)
        q=h['exchange'];t=h['t'];integral=np.r_[0.,np.cumsum(.5*(q[1:]+q[:-1])*np.diff(t))]
        residual=np.max(np.abs(h['axion_energy']-h['axion_energy'][0]-integral))/h['total'][0]
        results.append(dict(dt=dt,frequency_bound_absolute=float(abs(h['dt_frequency_bound'][-1]-expected_frequency)),state_error=float(np.max(np.abs(s-ref))),
            energy_drift=float(np.max(np.abs(h['total']-h['total'][0]))/h['total'][0]),
            exchange_residual=float(residual),exchange_amplitude=float(np.max(np.abs(q))),
            snapshot_energy_relative=float(abs(energy-h['total'][-1])/h['total'][0]),
            transfer_exact=bool(np.array_equal(state,np.loadtxt(directory/'native/initial.txt')))))
    translated,translated_final=execute(np.roll(state,7,axis=0),.0025,args.output/'translated')
    translation_error=max(float(np.max(np.abs(translated_final-np.roll(finals[-1],7,axis=0)))),
                          max(float(np.max(np.abs(translated[name]-histories[-1][name]))) for name in histories[-1].dtype.names))
    fine=results[-1]
    checks=dict(frequency=all(r['frequency_bound_absolute']<gates['frequency_bound_absolute'] for r in results),state=fine['state_error']<gates['finest_state_error'],
        second_order=all(results[j]['state_error']/results[j+1]['state_error']>gates['min_refinement_ratio'] for j in (0,1)),
        energy=fine['energy_drift']<gates['finest_energy_drift'],exchange=fine['exchange_residual']<gates['finest_exchange_residual'],
        snapshot_energy=all(r['snapshot_energy_relative']<gates['snapshot_energy_relative'] for r in results),
        transfer=all(r['transfer_exact'] for r in results),translation=translation_error<gates['translation_absolute'],
        nonzero_exchange=all(r['exchange_amplitude']>1e-5 for r in results))
    report=dict(kind=__doc__,scope='NDim=1 canonical periodic planar fields, energy per transverse area; not 3D stability',
        spatial_order=args.spatial_order,grid=n,box=box,duration=duration,time='Minkowski cosmic time',gates=gates,build=build,
        reference=dict(method='DOP853',rtol=2e-12,atol=2e-14,nfev=reference.nfev,spatial_operator=f'periodic 1D order-{args.spatial_order} stencil'),
        results=results,translation_absolute=translation_error,checks=checks,passed=all(checks.values()),
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('validate_line_backend.py','paper_reference.py','planar_stencils.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(passed=report['passed'],results=results,translation=translation_error,checks=checks),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
