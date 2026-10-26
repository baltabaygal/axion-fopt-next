"""Bounded native coupled-field distributed-slab pilot; not lifetime evidence."""
import argparse,json,hashlib,os,subprocess,time
from pathlib import Path
from dataclasses import replace
import numpy as np
from planar_distributed_initial import SlabSeed,prepare
from planar_diagnostics import wall_geometry

HERE=Path(__file__).resolve().parent

def integral(y,t):return np.r_[0.,np.cumsum(.5*(y[1:]+y[:-1])*np.diff(t))]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--executable',type=Path,required=True)
    ap.add_argument('--analyze-existing',action='store_true')
    ap.add_argument('--n',type=int,default=4096);ap.add_argument('--dt',type=float,default=.0125)
    args=ap.parse_args()
    if not 128<=args.n<=16384 or not np.isfinite(args.dt) or args.dt<=0 or abs(.25/args.dt-round(.25/args.dt))>1e-8:raise ValueError('Invalid resolution or history cadence')
    if not args.analyze_existing:args.output.mkdir(parents=True,exist_ok=False)
    elif not (args.output/'contract.json').exists() or (args.output/'report.json').exists():raise RuntimeError('Require unfinished saved pilot')
    exe=args.executable.resolve();build=json.loads((exe.parent/'build.json').read_text())
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Binary mismatch')
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale native build')
    gates=dict(solver_seconds=15.,energy=.002,exchange=.002,window_balance=.002,symmetry=1e-10,width_cells=4.,geometry=1e-12)
    contract=dict(gates=gates,n=args.n,box=512.,dt=args.dt,duration=320.,cadence=.25,
        delta_v=.0002,seed=SlabSeed().__dict__,window=dict(centre=256.,inner=144.,outer=160.,velocity=0.),
        hypothesis='Distributed fields reduce slab contraction relative to matched empty initial data',
        survival_claim=False,history_length_policy='2*(rising wall_x-256), conditional on reflection symmetry; both endpoints independently measured')
    if not args.analyze_existing:(args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    elif json.loads((args.output/'contract.json').read_text())!=contract:raise RuntimeError('Saved contract differs')
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    rows=[];hist={};spent=0.
    for label,ratio in [('empty',0.),('filled',1.)]:
        directory=args.output/label
        if not args.analyze_existing:directory.mkdir()
        x,s,m=prepare(args.n,replace(SlabSeed(),pressure_ratio=ratio))
        if not args.analyze_existing:
            np.savetxt(directory/'input.txt',s,fmt='%.17g')
            params=dict(N=args.n,lSide=512.,dt=args.dt,tMax=320.,expansion=False,evolver='PV2',PS_type=2,
                line_spatial_order=4,line_wall_diagnostics=True,line_window_centre=256.,line_window_inner=144.,line_window_outer=160.,line_window_velocity=0.,
                **{'lambda':8.},delta_v=.0002,mass_fv=.1,mass_kappa=float(np.log(100)),axion_fa=10.,output_every=round(.25/args.dt),
                initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
            (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
            remaining=15.-spent
            if remaining<=0:raise RuntimeError('Runtime cap exhausted')
            start=time.monotonic();run=subprocess.run([str(exe),'input='+str((directory/'input.in').resolve())],env=env,capture_output=True,text=True,timeout=remaining)
            spent+=time.monotonic()-start;(directory/'run.log').write_text(run.stdout+run.stderr)
            if run.returncode:raise RuntimeError('Native evolution failed; see log')
        h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True);f=np.loadtxt(directory/'native/final.txt');hist[label]=h
        right=wall_geometry(f,512.);left=512.-wall_geometry(f[::-1],512.)['position']
        sym=max(float(np.max(abs(f[:,[0,2]]-f[::-1,[0,2]]))),float(np.max(abs(f[:,[1,3]]+f[::-1,[1,3]]))))
        rows.append(dict(label=label,final_left=left,final_right=right['position'],final_length=right['position']-left,
            final_right_velocity=right['phase_velocity'],minimum_inferred_length=float(np.min(2*(h['wall_x']-256))),
            endpoint_symmetry=sym,minimum_width_cells=float(h['wall_width_cells'].min()),
            geometry_error=abs(right['position']-h['wall_x'][-1]),modes=m))
        print(label,'finished',spent,flush=True)
    h=hist['filled'];empty=hist['empty'];t=h['t'];e0=h['axion_energy'][0]
    matched=(h['phi_energy']-h['phi_energy'][0])-(empty['phi_energy']-empty['phi_energy'][0])
    closure=float(np.max(abs(matched+h['axion_energy']-e0))/e0)
    exchange=float(np.max(abs(h['axion_energy']-e0-integral(h['exchange'],t)))/e0)
    residual=h['window_axion_energy']-h['window_axion_energy'][0]+integral(h['window_axion_flux'],t)-integral(h['window_exchange'],t)
    window=float(np.max(abs(residual))/e0)
    checks=dict(energy=closure<gates['energy'],exchange=exchange<gates['exchange'],window=window<gates['window_balance'],
        symmetry=all(r['endpoint_symmetry']<gates['symmetry'] for r in rows),resolution=all(r['minimum_width_cells']>=4 for r in rows),
        geometry=all(r['geometry_error']<gates['geometry'] for r in rows),empty_exact=bool(np.all(empty['axion_energy']==0)))
    if args.analyze_existing:
        spent=float((args.output/'execution_stdout.txt').read_text().split('finished ')[-1].strip())
    result=dict(contract=contract,analysis_recovered_from_saved_fields=args.analyze_existing,rows=rows,checks=checks,numerical_gates_passed=all(checks.values()),
        reduced_contraction=bool(rows[1]['final_length']>rows[0]['final_length']),matched_energy_closure=closure,exchange_residual=exchange,
        fixed_window_balance_residual=window,final_axion_energy_ratio=float(h['axion_energy'][-1]/e0),
        fixed_window_net_outward_axion_flux_over_initial=float(integral(h['window_axion_flux'],t)[-1]/e0),
        total_solver_seconds=spent,survival_established=False,build=build,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('run_distributed_slab_pilot.py','planar_distributed_initial.py','planar_diagnostics.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('build','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
