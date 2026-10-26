"""Longer single bubble with frozen physics and explicit numerical controls."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from model import PARAMS,seed,read_state,run
from run_bubble import analyze
from validate_backend import verify_build
from build_backend import sha

CASES=[('main',128,48.,.04,1.),('coarse',96,48.,.04,1.),
       ('time-half',128,48.,.02,1.),('empty',128,48.,.04,0.),('larger-box',160,60.,.04,1.)]

def input_registry():
    return dict(physical=dict(c=1.,S_excursion_v=1.,canonical_fields=['S','A=f*theta'],
        potential='lambda*S^2*(1-S)^2/4 + DeltaV*(1-4*S^3+3*S^4) + f^2*M^2*S^2*(1-cos(A/f))',
        **PARAMS,scale_factor=1.,H=0.,friction=0.,temperature_model=None,plasma=False,
        backreaction=True,prescribed_wall_speed=None,stochastic_nucleation=False),
        initial=dict(profile='S=(1-tanh((r-R0)/w))/2',R0=8.,w=2.,centre='(L/2,L/2,L/2)',
            theta=1.,S_dot=0.,theta_dot=0.,random_seed=None,birth='one imposed seed at t=0'),
        numerical=dict(dimension=3,boundaries='periodic on each axis',time='Minkowski cosmic time',
            t_initial=0.,t_final=18.,cases=[dict(name=a,N=b,L=c,dx=c/b,dt=d,theta=e) for a,b,c,d,e in CASES],
            evolver='CosmoLattice PV2',laplacian='nearest-neighbor second order',
            energy_gradient='compatible forward difference',precision='float64',processes=1,threads=1,
            history_cadence=.08,slice_cadence=.4,full_states='initial and final',subprocess_timeout_seconds=180,
            withGWs=False,fixedBackground=True,omegaEoS=0,PS_type=1,
            kIR='2*pi/L; metadata, no spectrum estimator used',
            native_files='input.bin, input.in -> native/history.csv, frames and initial/final states'))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--validation',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--analyze-only',action='store_true');args=p.parse_args()
    out=args.output.resolve();build=verify_build(args.executable)
    val=json.loads(args.validation.read_text())
    if not val['passed'] or val['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Matching validated executable required')
    if not args.analyze_only:
        out.mkdir(parents=True,exist_ok=False);shutil.copy2(Path(__file__).with_name('LONG_EXPERIMENT.md'),out/'contract.md')
        (out/'parameters.json').write_text(json.dumps(input_registry(),indent=2)+'\n')
        (out/'build.json').write_text(json.dumps(build,indent=2)+'\n');sources=out/'sources';sources.mkdir()
        for f in Path(__file__).parent.iterdir():
            if f.suffix in ('.py','.cpp','.h','.md'):shutil.copy2(f,sources/f.name)
        timings=[]
        for name,n,box,dt,theta in CASES:
            sec=run(args.executable,out/name,seed(n,box,theta),box,dt,18.,timeout=180)
            item=dict(case=name,seconds=sec);timings.append(item)
            (out/'timings.json').write_text(json.dumps(timings,indent=2)+'\n');print(json.dumps(item),flush=True)
    results={name:analyze(out/name,n,box,dt,theta) for name,n,box,dt,theta in CASES}
    main=results['main'];coarse=results['coarse'];half=results['time-half'];large=results['larger-box']
    _,y=read_state(out/'main/native/final.bin');_,yh=read_state(out/'time-half/native/final.bin')
    scale=np.array([1.,PARAMS['fa'],1.,PARAMS['fa']]);time_rms=float(np.sqrt(np.mean(((y-yh)/scale)**2)));del yh
    _,yl=read_state(out/'larger-box/native/final.bin');aligned=yl[16:-16,16:-16,16:-16]
    x=(np.arange(128)+.5)*.375-24;mask=x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2<20**2
    boundary_rms=float(np.sqrt(np.mean(((y[mask]-aligned[mask])/scale)**2)))
    relative=lambda a,b:float(abs(a-b)/max(abs(b),1e-30)) if a is not None and b is not None else None
    metrics=dict(time_field_rms=time_rms,time_radius=relative(main['final_radius'],half['final_radius']),
        time_axion_energy=relative(main['final_axion_energy'],half['final_axion_energy']),
        spatial_radius=relative(coarse['final_radius'],main['final_radius']),
        spatial_axion_energy=relative(coarse['final_axion_energy'],main['final_axion_energy']),
        boundary_field_rms=boundary_rms,boundary_radius=relative(main['final_radius'],large['final_radius']),
        boundary_global_axion_energy=relative(main['final_axion_energy'],large['final_axion_energy']))
    wall_checks={}
    for name in ('main','empty'):
        r=results[name];dx=r['box']/r['n'];fail=[];widths=[]
        for c in r['wall_trajectory']:
            if c['width'] is not None:widths.append(c['width']/dx)
            if c['width'] is None or c['width']/dx<5 or any(v!=1 for v in c['crossing_counts'].values()):fail.append(c['t'])
        wall_checks[name]=dict(passed=not fail,min_cells=min(widths) if widths else None,first_failure_time=fail[0] if fail else None)
    checks=dict(reference=val['passed'],initial_transfer=all(r['initial_transfer_exact'] for r in results.values()),
        energy_drift=all(r['energy_drift_excess']<.005 for r in results.values()),
        snapshot_energy=all(r['snapshot_energy_relative']<1e-11 for r in results.values()),
        exchange=all(r['exchange_residual']<.01 for r in results.values() if r['theta_initial']!=0),
        empty_exact=results['empty']['empty_exact'],time_fields=time_rms<.01,
        time_radius=metrics['time_radius'] is not None and metrics['time_radius']<.005,
        time_energy=metrics['time_axion_energy']<.01,
        spatial_radius=metrics['spatial_radius'] is not None and metrics['spatial_radius']<.05,
        spatial_energy=metrics['spatial_axion_energy']<.1,
        resolved_walls=all(r['passed'] for r in wall_checks.values()),
        boundary_fields=boundary_rms<.001,boundary_radius=metrics['boundary_radius'] is not None and metrics['boundary_radius']<.005)
    report=dict(passed=bool(all(checks.values())),checks=checks,metrics=metrics,wall_checks=wall_checks,results=results,
        build=build,parameters=input_registry(),validation_sha256=sha(args.validation),
        source_sha256={str(f.resolve()):sha(f) for f in Path(__file__).parent.iterdir() if f.suffix in ('.py','.cpp','.h','.md')},
        scope='Finite-time seeded 3D coupled fields, not a continuum limit, terminal speed or cosmological abundance')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k in ('passed','checks','metrics','wall_checks')},indent=2),flush=True)
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
