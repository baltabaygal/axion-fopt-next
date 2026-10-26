"""Periodic small-box completion pilot; fresh simulations and global diagnostics."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from model import PARAMS,seed,read_state,run,terms
from run_bubble import contour
from validate_backend import verify_build
from build_backend import sha

CASES=[('main',96,.02),('coarse',64,.02),('time-half',96,.01)]
BOX=24.;END=32.

def phase_values(s):
    return dict(S_min=float(s.min()),S_max=float(s.max()),S_mean=float(s.mean()),S_std=float(s.std()),
        F_half=float(np.mean(s<.5)),F_barrier=float(np.mean(s<2/(4+24*.06))),near_true_fraction=float(np.mean(abs(s-1)<.1)))

def intervals(t,mask):
    blocks=[];begin=None
    for j,yes in enumerate(mask):
        if yes and begin is None:begin=j
        if begin is not None and (not yes or j==len(mask)-1):
            stop=j-1 if not yes else j
            blocks.append(dict(start=float(t[begin]),end=float(t[stop]),duration=float(t[stop]-t[begin])));begin=None
    return blocks

def analyze(directory,n,dt,box=BOX,end=END):
    native=directory/'native';h=np.genfromtxt(native/'history.csv',delimiter=',',names=True);dx=box/n
    _,initial=read_state(native/'initial.bin');_,final=read_state(native/'final.bin')
    snapshots=[]
    for y,j in [(initial,0),(final,-1)]:
        phase=phase_values(y[...,0]);energies={k:float(v.sum()*dx**3) for k,v in terms(y,dx).items()}
        snapshots.append(dict(phase=phase,phase_moment_error=max(abs(v-h[k][j]) for k,v in phase.items() if k.startswith('S_')),
            fractions_exact=all(v==h[k][j] for k,v in phase.items() if not k.startswith('S_')),
            energy_error=max(abs(v-h[k][j]) for k,v in energies.items())/abs(h['total'][0])))
    full=h['F_half']==0;blocks=intervals(h['t'],full);sustained=[v for v in blocks if v['duration']>=2]
    first=float(h['t'][np.flatnonzero(full)[0]]) if full.any() else None
    recross=bool(first is not None and np.any(h['F_half'][h['t']>first]>0))
    q=np.r_[0,np.cumsum((h['Q_A'][1:]+h['Q_A'][:-1])/2*np.diff(h['t']))]
    initial_excess=float(h['total'][0]-.06*box**3)
    drift=np.max(abs(h['total']-h['total'][0]))
    wall=[];stop=None
    for path in sorted(native.glob('frame_*.bin')):
        meta,plane=read_state(path);c=contour(plane[:,n//2-1:n//2+1,0].mean(axis=1),box,dx*dx/2)
        if c['width'] is None or c['r10']>=box/2:
            stop=meta['t'];break
        wall.append(dict(t=meta['t'],width_cells=c['width']/dx,radius=c['radius']))
    checks=dict(initial_exact=bool(np.array_equal(initial,seed(n,box,1.))),finite=bool(np.isfinite(final).all()),
        independent_energy=all(v['energy_error']<1e-11 for v in snapshots),
        independent_phase=all(v['fractions_exact'] and v['phase_moment_error']<1e-12 for v in snapshots),
        energy_drift=float(drift/abs(initial_excess))<.005,
        exchange=float(np.max(abs(h['E_A']-h['E_A'][0]-q))/h['E_A'][0])<.01,
        pre_collision_wall=bool(wall) and min(v['width_cells'] for v in wall)>=5)
    return dict(n=n,box=box,dt=dt,end=end,checks=checks,snapshots=snapshots,
        first_full_conversion=first,first_sustained_full_interval=sustained[0] if sustained else None,
        full_intervals=blocks,recrossing_after_first=recross,final_phase=phase_values(final[...,0]),
        initial_excess=initial_excess,energy_drift_excess=float(drift/abs(initial_excess)),energy_drift_total=float(drift/abs(h['total'][0])),
        exchange_residual=float(np.max(abs(h['E_A']-h['E_A'][0]-q))/h['E_A'][0]),
        final_axion_energy=float(h['E_A'][-1]),pre_collision_wall=wall,radial_diagnostic_stops_at=stop)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--validation',type=Path,required=True);p.add_argument('--previous-validation',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output.resolve()
    build=verify_build(args.executable);val=json.loads((args.validation/'report.json').read_text())
    if not val['passed'] or val['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Validated matching backend required')
    out.mkdir(parents=True,exist_ok=False);sources=out/'sources';sources.mkdir()
    shutil.copy2(Path(__file__).with_name('SMALL_BOX_EXPERIMENT.md'),out/'contract.md')
    (out/'build.json').write_text(json.dumps(build,indent=2)+'\n')
    for f in Path(__file__).parent.iterdir():
        if f.suffix in ('.py','.h','.cpp','.md'):shutil.copy2(f,sources/f.name)
    regression=[]
    for case in args.previous_validation.iterdir():
        if case.is_dir() and (case/'native/final.bin').exists():
            _,old=read_state(case/'native/final.bin');_,new=read_state(args.validation/case.name/'native/final.bin')
            regression.append(bool(np.array_equal(old,new)))
    if not regression or not all(regression):raise RuntimeError('Passive diagnostics changed reference fields')
    settings=dict(physical=PARAMS,box=BOX,end=END,initial=dict(R0=8.,width=2.,theta=1.,S_dot=0.,theta_dot=0.),
        cases=[dict(name=a,n=b,dx=BOX/b,dt=c) for a,b,c in CASES],cadence=.04,slice_cadence=.2,
        boundary='periodic; image interaction is part of this experiment',completion='all 3D cells S>=.5, measured every .04',
        no_expansion=True,no_friction=True,backreaction=True)
    (out/'parameters.json').write_text(json.dumps(settings,indent=2)+'\n')
    results={};timing=[]
    for name,n,dt in CASES:
        sec=run(args.executable,out/name,seed(n,BOX,1.),BOX,dt,END,cadence=.04,frames=.2,timeout=240)
        r=analyze(out/name,n,dt);results[name]=r
        (out/name/'analysis.json').write_text(json.dumps(r,indent=2)+'\n')
        timing.append(dict(case=name,seconds=sec));(out/'timings.json').write_text(json.dumps(timing,indent=2)+'\n')
        print(json.dumps(dict(case=name,seconds=sec,first_full=r['first_full_conversion'],sustained=r['first_sustained_full_interval'],final_phase=r['final_phase'],checks=r['checks'])),flush=True)
    base=results['main'];coarse=results['coarse'];half=results['time-half']
    def compare(r):
        a,b=r['first_full_conversion'],base['first_full_conversion']
        return dict(energy_relative=abs(r['final_axion_energy']/base['final_axion_energy']-1),
            false_fraction=abs(r['final_phase']['F_half']-base['final_phase']['F_half']),
            completion_time_difference=abs(a-b) if a is not None and b is not None else None)
    spatial=compare(coarse);temporal=compare(half)
    checks=dict(reference=val['passed'],exact_reference_regression=all(regression),
        local_checks=all(all(r['checks'].values()) for r in results.values()),
        spatial_energy=spatial['energy_relative']<.05,spatial_phase=spatial['false_fraction']<.01,
        spatial_completion=spatial['completion_time_difference'] is not None and spatial['completion_time_difference']<=.4,
        temporal_energy=temporal['energy_relative']<.01,temporal_phase=temporal['false_fraction']<.005,
        temporal_completion=temporal['completion_time_difference'] is not None and temporal['completion_time_difference']<=.08)
    report=dict(passed=bool(all(checks.values())),checks=checks,results=results,spatial=spatial,temporal=temporal,
        goal_reached=base['first_sustained_full_interval'] is not None,parameters=settings,build=build,
        validation_sha256=sha(args.validation/'report.json'),source_sha256={str(f.resolve()):sha(f) for f in sources.iterdir()},
        caveat='Operational true phase does not mean uniform static S=1 or constant axion mass. Periodic-image collisions are intentional.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k in ('passed','checks','goal_reached','spatial','temporal')},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
