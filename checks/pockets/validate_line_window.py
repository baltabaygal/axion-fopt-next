"""Tiny compiled passive-window agreement and integrated continuity gates."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from paper_reference import Parameters
from planar_stencils import window_budget

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--spatial-order',type=int,choices=(2,4),default=2)
    parser.add_argument('--baseline',type=Path,required=True,help='Saved formation-line-validation-01 directory')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale line build')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    if args.spatial_order not in build.get('supports_line_spatial_orders',[2]):raise RuntimeError('Unsupported compiled spatial order')
    previous=json.loads((args.baseline/'report.json').read_text())
    if previous.get('spatial_order',2)!=args.spatial_order:raise RuntimeError('Baseline spatial order mismatch')
    gates=dict(reference_relative=1e-12,finest_continuity=1e-5,minimum_refinement_ratio=3.,translation_absolute=1e-12)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    p=Parameters();n=32;box=8.;dx=box/n;x=(np.arange(n)+.5)*dx
    source=args.baseline/'dt-0.01';state=np.loadtxt(source/'native/initial.txt')
    def execute(s,dt,directory,shift=0):
        directory.mkdir();np.savetxt(directory/'input.txt',s,fmt='%.17g')
        params=dict(line_spatial_order=args.spatial_order,N=n,lSide=box,dt=dt,tMax=.2,expansion=False,evolver='PV2',PS_type=2,
            **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,output_every=1,
            line_window_inner=.5,line_window_outer=1.5,line_window_centre=(7.7+shift*dx)%box,line_window_velocity=.2,
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
        run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],env=environment,capture_output=True,text=True,timeout=10)
        (directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        return np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True),np.loadtxt(directory/'native/final.txt')
    rows=[];histories=[];finals=[]
    for dt in (.01,.005,.0025):
        directory=args.output/f'dt-{dt}';h,final=execute(state,dt,directory);histories.append(h);finals.append(final)
        max_agreement=0.
        for j,s,t in ((0,state,0.),(-1,final,.2)):
            d=(x-7.7-.2*t+box/2)%box-box/2;z=np.clip((abs(d)-.5),0,1)
            w=1-z**3*(10+z*(-15+6*z));wd=.2*np.sign(d)*30*z*z*(1-z)**2
            ref=window_budget(s,box,w,wd,p,args.spatial_order)
            (directory/f'reference-{j}.json').write_text(json.dumps(ref,indent=2)+'\n')
            max_agreement=max(max_agreement,max(abs(v-h['window_'+k][j]) for k,v in ref.items())/h['window_energy'][0])
        t=h['t'];integrate=lambda y:np.r_[0.,np.cumsum(.5*(y[1:]+y[:-1])*np.diff(t))]
        q=integrate(h['window_exchange']);fp=integrate(h['window_phi_flux']);fa=integrate(h['window_axion_flux'])
        residuals=dict(total=h['window_energy']-h['window_energy'][0]+fp+fa,
            phi=h['window_phi_energy']-h['window_phi_energy'][0]+fp+q,
            axion=h['window_axion_energy']-h['window_axion_energy'][0]+fa-q)
        residual=max(float(np.max(abs(v))/h['window_energy'][0]) for v in residuals.values())
        previous=args.baseline/f'dt-{dt}'/'native'
        old=np.genfromtxt(previous/'history.csv',delimiter=',',names=True)
        unchanged=bool((directory/'native/final.txt').read_bytes()==(previous/'final.txt').read_bytes() and
            all(np.array_equal(h[name],old[name]) for name in old.dtype.names))
        rows.append(dict(dt=dt,reference_relative=float(max_agreement),continuity_relative=residual,
            passive_fields_and_old_history_exact=unchanged,exchange_amplitude=float(np.max(abs(h['window_exchange']))),
            advection_amplitude=float(np.max(abs(h['window_phi_advection'])))))
    shifted,shifted_final=execute(np.roll(state,7,axis=0),.0025,args.output/'shifted',7)
    translation=max(float(np.max(abs(shifted_final-np.roll(finals[-1],7,axis=0)))),
        max(float(np.max(abs(shifted[name]-histories[-1][name]))) for name in shifted.dtype.names))
    checks=dict(reference=all(r['reference_relative']<gates['reference_relative'] for r in rows),
        integrated_continuity=rows[-1]['continuity_relative']<gates['finest_continuity'],
        convergence=all(rows[j]['continuity_relative']/rows[j+1]['continuity_relative']>gates['minimum_refinement_ratio'] for j in (0,1)),
        passive=all(r['passive_fields_and_old_history_exact'] for r in rows),
        translation=translation<gates['translation_absolute'],
        nontrivial=all(r['exchange_amplitude']>1e-6 and r['advection_amplitude']>1e-4 for r in rows))
    report=dict(kind=__doc__,spatial_order=args.spatial_order,grid=n,box=box,duration=.2,time='Minkowski cosmic time',
        window=dict(inner=.5,outer=1.5,centre=7.7,velocity=.2,description='Periodic compact translating diagnostic, no force modification'),
        gates=gates,rows=rows,translation_absolute=translation,checks=checks,passed=all(checks.values()),build=build,
        baseline=str(args.baseline.resolve()),source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('validate_line_window.py','planar_stencils.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],rows=rows,translation=translation,checks=checks),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
