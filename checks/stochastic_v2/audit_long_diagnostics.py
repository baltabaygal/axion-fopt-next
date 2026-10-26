"""Repair only endpoint accounting via compensated native diagnostic replay.

Preserves the original failed report and all original evolution histories.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from model import read_state,run,terms
from build_backend import sha
from validate_backend import verify_build

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--old-validation',type=Path,required=True);p.add_argument('--new-validation',type=Path,required=True)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    build=verify_build(args.executable);original=json.loads((args.source/'report.json').read_text())
    newval=json.loads((args.new_validation/'report.json').read_text())
    if not newval['passed'] or newval['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Validated matching executable required')
    # This is a predeclared equality audit, not a replacement physical reference.
    (out/'contract.json').write_text(json.dumps(dict(snapshot_energy_gate=1e-11,exact_state_regression=True,
        scope='Replay initial diagnostics of saved endpoints; no long trajectory is replaced'),indent=2)+'\n')
    regression=[]
    for directory in args.old_validation.iterdir():
        if not directory.is_dir() or not (directory/'native/final.bin').exists():continue
        _,old=read_state(directory/'native/final.bin');_,new=read_state(args.new_validation/directory.name/'native/final.bin')
        regression.append(dict(case=directory.name,exact=bool(np.array_equal(old,new))))
    # Check that only the accumulator and its include changed in the C++ driver.
    old_source=(args.source/'sources/evolve.cpp').read_text()
    expected=old_source.replace('#include "../pockets/initial_state.h"','#include "../pockets/initial_state.h"\n#include "../pockets/accurate_sum.h"').replace(
        'long double es=0,ea=0,exchange=0,ks=0,gs=0,us=0,ka=0,ga=0,va=0;','AccurateSum es,ea,exchange,ks,gs,us,ka,ga,va;')
    diagnostic_only=expected==Path(__file__).with_name('evolve.cpp').read_text()
    snapshots=[]
    for name,r in original['results'].items():
        meta,y=read_state(args.source/name/'native/final.bin')
        sec=run(args.executable,out/name,y,meta['box'],.02,.02,cadence=.02,frames=.02,timeout=60)
        _,transferred=read_state(out/name/'native/initial.bin')
        h=np.genfromtxt(out/name/'native/history.csv',delimiter=',',names=True)
        energies={k:float(v.sum()*(meta['box']/meta['n'])**3) for k,v in terms(y,meta['box']/meta['n']).items()}
        error=max(abs(v-h[k][0]) for k,v in energies.items())/abs(r['initial_energy'])
        snapshots.append(dict(case=name,initial_transfer_exact=bool(np.array_equal(y,transferred)),
            corrected_native_vs_independent_relative=float(error),original_relative=r['snapshot_energy_relative'],seconds=sec))
    checks=dict(original_evolution_checks=all(v for k,v in original['checks'].items() if k!='snapshot_energy'),
        new_reference=newval['passed'],diagnostic_only_source_change=diagnostic_only,
        small_evolution_exact_regression=bool(regression) and all(v['exact'] for v in regression),
        corrected_snapshot_energy=all(v['initial_transfer_exact'] and v['corrected_native_vs_independent_relative']<1e-11 for v in snapshots))
    report=dict(passed=bool(all(checks.values())),checks=checks,original_checks=original['checks'],snapshots=snapshots,
        exact_field_regression=regression,original_report_sha256=sha(args.source/'report.json'),build=build,
        script_sha256=sha(__file__),scope=__doc__)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='build'},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
