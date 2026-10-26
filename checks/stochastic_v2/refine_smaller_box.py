"""Validate the observed transient all-cell interval, preserving original failures."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from model import seed,run
from run_small_box import analyze
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();build=verify_build(args.executable);out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(Path(__file__).with_name('SMALL_BOX_REFINEMENT.md'),out/'contract.md')
    base=json.loads((args.baseline/'report.json').read_text())['results'];bh=np.genfromtxt(args.baseline/'main/native/history.csv',delimiter=',',names=True)
    if not base['full_intervals']:raise RuntimeError('No baseline transient interval')
    interval=base['full_intervals'][0];base_e=float(np.interp(30,bh['t'],bh['E_A']))
    results={};timings=[]
    for name,n,dt,shift,energy_gate in [('spatial',96,.01,.4,.05),('temporal',80,.005,.08,.01)]:
        seconds=run(args.executable,out/name,seed(n,20.,1.),20.,dt,30.,cadence=.04,frames=.2,timeout=180)
        r=analyze(out/name,n,dt,box=20.,end=30.);blocks=r['full_intervals']
        closest=min(blocks,key=lambda v:abs(v['start']-interval['start'])) if blocks else None
        checks={k:v for k,v in r['checks'].items() if k!='pre_collision_wall'}
        checks.update(transient_interval=closest is not None,
            interval_start=closest is not None and abs(closest['start']-interval['start'])<=shift,
            interval_end=closest is not None and abs(closest['end']-interval['end'])<=shift,
            energy_comparison=abs(r['final_axion_energy']/base_e-1)<energy_gate)
        results[name]=dict(analysis=r,checks=checks,passed=bool(all(checks.values())),energy_difference=abs(r['final_axion_energy']/base_e-1),
            radius_observable='not applicable: periodic interaction from initial tails')
        (out/(name+'_report.json')).write_text(json.dumps(results[name],indent=2)+'\n')
        timings.append(dict(case=name,seconds=seconds));(out/'timings.json').write_text(json.dumps(timings,indent=2)+'\n')
        print(json.dumps(dict(case=name,seconds=seconds,checks=checks,full_intervals=blocks)),flush=True)
    report=dict(passed=bool(all(r['passed'] for r in results.values())),results=results,baseline_interval=interval,
        original_sustained_goal_reached=base['first_sustained_full_interval'] is not None,
        build=build,baseline_report_sha256=sha(args.baseline/'report.json'),
        source_sha256={str(f.resolve()):sha(f) for f in [Path(__file__),Path(__file__).with_name('model.py'),Path(__file__).with_name('run_small_box.py')]},
        scope='Transient all-cell phase-threshold crossing only; original failed reports are preserved')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(passed=report['passed'],original_sustained_goal_reached=report['original_sustained_goal_reached'])))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
