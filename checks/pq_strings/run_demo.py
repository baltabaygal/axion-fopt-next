"""Run deterministic/refinement and small cooled-network PQ experiments."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import time
import numpy as np
import scipy
from model import (ROOT, EXECUTABLE, config, run, read_state, coefficients,
                   deterministic, gaussian, cool, defects, energies)

# Stored before any run; failed checks remain in the report.
CONTRACT = dict(
    model='arXiv:2609.20801v1 Eq.3; two real canonical fields; f_a=1',
    background='radiation a=tau/20, pi=a^2*F_prime',
    lambda_=0.08, tau_origin=20., tau_star=27., n_tau=8.,
    preparation='random Gaussian ell=2.5, seed=0, E[|Phi|^2]=1; static chi=0 diffusion time 2',
    simplification='no core-growth stage; short finite-box evolution; no scaling or abundance inference',
    gates=dict(minimum_core_cells=4., minimum_wall_cells=4., maximum_ma_ms=.16,
               exact_initial_transfer=True, snapshot_energy_error=2e-12,
               balance_relative=.003, temporal_field_rms=.001,
               spatial_field_rms=.05, cooling_energy_increase=1e-12))


def summarize(directory, cfg, initial):
    directory=Path(directory);native=directory/'native'
    h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
    _,actual=read_state(native/'initial.bin')
    gates=dict(initial_transfer=bool(np.array_equal(actual,initial)),finite_history=bool(np.isfinite(h.tolist()).all()))
    snapshots=[];max_energy_error=0.
    for path in sorted(native.glob('frame_*.bin')):
        meta,s=read_state(path);diag,_=defects(s,meta['box'],meta['tau'])
        e=energies(s,meta['box'],meta['a'],cfg['lambda'],meta['chi'])
        matching=np.where(abs(h['u']-meta['u'])<1e-10)[0]
        if len(matching):
            max_energy_error=max(max_energy_error,max(abs(e[k]-h[k][matching[0]]) for k in e))
        ms=np.sqrt(2*cfg['lambda']);ma=np.sqrt(meta['chi']);dx=meta['box']/meta['n']
        row=dict(path=path.name,**meta,**diag,**e,
                 core_cells=1/(meta['a']*ms*dx),
                 wall_cells=None if ma==0 else 1/(meta['a']*ma*dx),ma_ms=ma/ms,
                 r_min=float(np.linalg.norm(s[...,:2],axis=-1).min()),
                 r_mean=float(np.linalg.norm(s[...,:2],axis=-1).mean()))
        snapshots.append(row)
    meta,final=read_state(native/'final.bin')
    diag,_=defects(final,meta['box'],meta['tau'])
    core_min=1/(meta['a']*np.sqrt(2*cfg['lambda'])*(meta['box']/meta['n']))
    wall_min=None if meta['chi']==0 else 1/(meta['a']*np.sqrt(meta['chi'])*(meta['box']/meta['n']))
    mass_ratio=np.sqrt(meta['chi']/(2*cfg['lambda']))
    work=np.trapezoid(h['source']+h['dilution'],h['u'])
    residual=abs(h['total'][-1]-h['total'][0]-work)/max(h['total'][0],h['total'][-1],1e-12)
    gates.update(core_resolution=core_min>=CONTRACT['gates']['minimum_core_cells'],
                 wall_resolution=wall_min is None or wall_min>=CONTRACT['gates']['minimum_wall_cells'],
                 mass_hierarchy=mass_ratio<=CONTRACT['gates']['maximum_ma_ms'],
                 energy_transfer=max_energy_error<CONTRACT['gates']['snapshot_energy_error'],
                 energy_balance=residual<CONTRACT['gates']['balance_relative'])
    report=dict(parameters=cfg,gates={k:bool(v) for k,v in gates.items()},
                passed=all(gates.values()),minimum_core_cells=core_min,minimum_wall_cells=wall_min,
                maximum_ma_ms=mass_ratio,energy_balance_relative=float(residual),
                snapshot_energy_error=max_energy_error,initial_defects=defects(initial,cfg['lSide'],cfg['tau_origin'])[0],
                final_defects=diag,final_energy=float(h['total'][-1]),snapshots=snapshots)
    (directory/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return report,final


def execute(name,initial,cfg,root,executable):
    started=time.perf_counter()
    actual=run(initial,cfg,root/name,executable)
    report,final=summarize(root/name,actual,initial)
    report['elapsed_seconds']=time.perf_counter()-started
    (root/name/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(name, 'passed', report['passed'], 'plaquettes',report['initial_defects']['plaquettes'],
          '->',report['final_defects']['plaquettes'],'seconds',round(report['elapsed_seconds'],2),flush=True)
    return report,final


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=EXECUTABLE)
    p.add_argument('--suite',choices=['pilot','refinement','network','matched','collapse'],default='pilot')
    args=p.parse_args();root=args.output;root.mkdir(parents=True,exist_ok=False)
    (root/'contract.json').write_text(json.dumps(CONTRACT,indent=2)+'\n')
    provenance=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
        sources={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in sorted((ROOT/'checks/pq_strings').glob('*')) if path.is_file()},
        build=json.loads((args.executable.parent/'adapter/build-info.json').read_text()))
    (root/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    shutil.copytree(ROOT/'checks/pq_strings',root/'source',ignore=shutil.ignore_patterns('__pycache__','results'))
    reports={};comparisons={};common=dict(lam=.08,origin=20.,star=27.,exponent=8.,every=5,frames=25)
    if args.suite in ('pilot','refinement'):
        n=32 if args.suite=='pilot' else 64
        d=deterministic(n,12.,core=2.5)
        for label,chi in [('growing',None),('massless',0.)]:
            cfg=config(n,12.,.02,6.,chi_star=chi,**common)
            reports['deterministic_'+label],final=execute('deterministic_'+label,d,cfg,root,args.executable)
            if label=='growing':coarse=final
        cfg=config(n,12.,.01,6.,**dict(common,every=10,frames=50))
        reports['deterministic_fine_dt'],fine=execute('deterministic_fine_dt',d,cfg,root,args.executable)
        rms=float(np.sqrt(np.mean((coarse[...,:2]-fine[...,:2])**2)))
        comparisons['temporal_field_rms']=rms
        comparisons['temporal_passed']=rms<CONTRACT['gates']['temporal_field_rms']
        if args.suite=='pilot':
            raw=gaussian(32,12.,2.5,0);s,cooling=cool(raw,12.,.08,duration=2.)
            np.save(root/'random_raw.npy',raw);np.save(root/'random_cooled.npy',s)
            np.savetxt(root/'cooling.csv',cooling,delimiter=',',header='diffusion_time,total',comments='')
            comparisons['cooling_passed']=bool(np.max(np.diff(cooling[:,1]))<CONTRACT['gates']['cooling_energy_increase'])
            for label,chi in [('growing',None),('massless',0.)]:
                cfg=config(32,12.,.02,6.,chi_star=chi,**common)
                reports['network_'+label],_=execute('network_'+label,s,cfg,root,args.executable)
    elif args.suite == 'collapse':
        comparisons['initial_definition']='Four neutral analytic straight strings; X=cos(2pi*x/L)+0.5, Y=sin(2pi*y/L), tanh modulus core=2.5, L=12'
        finals=[]
        for n in (32,64):
            initial=deterministic(n,12.,core=2.5,bias=.5)
            cfg=config(n,12.,.02,6.,**common)
            reports[f'pair_N{n}'],final=execute(f'pair_N{n}',initial,cfg,root,args.executable)
            finals.append(final)
        rms=float(np.sqrt(np.mean((finals[0][...,:2]-finals[1].reshape(32,2,32,2,32,2,4).mean(axis=(1,3,5))[...,:2])**2)))
        comparisons.update(spatial_field_rms=rms,spatial_passed=rms<CONTRACT['gates']['spatial_field_rms'])
        initial=deterministic(32,12.,core=2.5,bias=.5)
        cfg=config(32,12.,.02,6.,chi_star=0.,**common)
        reports['pair_massless'],_=execute('pair_massless',initial,cfg,root,args.executable)
    elif args.suite == 'matched':
        # Freeze one fine-grid preparation; coarse data are its cell-volume averages.
        raw=gaussian(64,12.,2.5,0);fine_initial,cooling=cool(raw,12.,.08,duration=2.)
        coarse_initial=fine_initial.reshape(32,2,32,2,32,2,4).mean(axis=(1,3,5))
        np.save(root/'random_raw.npy',raw);np.save(root/'random_cooled.npy',fine_initial)
        np.savetxt(root/'cooling.csv',cooling,delimiter=',',header='diffusion_time,total',comments='')
        comparisons['initial_matching']='N32 cell averages of frozen cooled N64 initial state, L=12'
        results=[]
        for n,initial in [(32,coarse_initial),(64,fine_initial)]:
            cfg=config(n,12.,.02,6.,**common)
            reports[f'network_N{n}'],final=execute(f'network_N{n}',initial,cfg,root,args.executable)
            results.append(final)
        block=results[1].reshape(32,2,32,2,32,2,4).mean(axis=(1,3,5))
        rms=float(np.sqrt(np.mean((results[0][...,:2]-block[...,:2])**2)))
        comparisons.update(spatial_field_rms=rms,spatial_passed=rms<CONTRACT['gates']['spatial_field_rms'])
    else:
        # Larger volume remains 64^3: dx=.375, same as the successful 32^3 pilot.
        raw=gaussian(64,24.,2.5,0);s,cooling=cool(raw,24.,.08,duration=2.)
        np.save(root/'random_raw.npy',raw);np.save(root/'random_cooled.npy',s)
        np.savetxt(root/'cooling.csv',cooling,delimiter=',',header='diffusion_time,total',comments='')
        comparisons['cooling_passed']=bool(np.max(np.diff(cooling[:,1]))<CONTRACT['gates']['cooling_energy_increase'])
        for label,chi in [('growing',None),('massless',0.)]:
            cfg=config(64,24.,.02,12.,chi_star=chi,**dict(common,frames=50))
            reports['network_'+label],final=execute('network_'+label,s,cfg,root,args.executable)
            if label=='growing':coarse=final
        cfg=config(64,24.,.01,12.,**dict(common,every=10,frames=100))
        reports['network_fine_dt'],fine=execute('network_fine_dt',s,cfg,root,args.executable)
        rms=float(np.sqrt(np.mean((coarse[...,:2]-fine[...,:2])**2)))
        comparisons.update(temporal_field_rms=rms,temporal_passed=rms<CONTRACT['gates']['temporal_field_rms'])
    passed=all(r['passed'] for r in reports.values()) and all(bool(v) for k,v in comparisons.items() if k.endswith('_passed'))
    result=dict(suite=args.suite,contract=CONTRACT,provenance=provenance,runs=reports,comparisons=comparisons,passed=passed)
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Suite passed:',passed,flush=True)
    raise SystemExit(0 if passed else 1)
if __name__=='__main__':main()
