"""Fixed-state field refinement; predeclared in BALANCED_REFINEMENT.md."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from model import (EXECUTABLE, config, eigen_state, history, radial_modes,
                   read_state, run)
from run_prototype import wall_cells
from model import read_slice
from wall_diagnostics import local_widths

RUNS = dict(time=dict(n=128, box=64., dt=.00625),
            space=dict(n=192, box=64., dt=.0125),
            box=dict(n=192, box=96., dt=.0125))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(name, root, baseline, modes):
    spec = RUNS[name];n, box, dt = (spec[k] for k in ('n', 'box', 'dt'))
    initial = eigen_state(n, box, 20., 2., modes, baseline['amplitude'], 0, 30.)
    cfg = config(n, box, dt, 80., lam=2., dv=.01, mass_tv=2.,
                 every=round(.5/dt), slices=round(4/dt))
    error = None
    try:
        run(initial, cfg, root/name)
    except RuntimeError as exc:
        error = str(exc)
    h = history(root/name)
    _, transferred = read_state(root/name/'native/initial.bin')
    e0 = sum(float(h[k][0]) for k in ('axion_kinetic','axion_gradient','coupling'))
    drift = float(np.max(np.abs(h['total']/h['total'][0]-1)))
    original, diagnostic = [], []
    for path in sorted((root/name/'native').glob('slice_*.bin')):
        meta, s = read_slice(path)
        w = wall_cells(s, box/n)
        original.append(dict(time=meta['u'], width=float(w) if np.isfinite(w) else None))
        diagnostic.append(dict(time=meta['u'], **local_widths(s[...,0], box/n)))
    gates = dict(transfer=bool(np.array_equal(initial, transferred)),
                 completed=error is None and bool(np.isclose(h['u'][-1],80.)),
                 energy=drift < 5e-4)
    report = dict(spec=spec, gates=gates, solver_error=error, initial_axion_energy=e0,
                  energy_drift=drift, final_radius=float(h['fv_radius'][-1]),
                  final_far_storage=float(h['axion_far'][-1]/baseline['target_axion_energy']),
                  original_wall=original, slice_normals=diagnostic,
                  slice_diagnostic_pass=all(d['diagnostic_pass'] for d in diagnostic))
    (root/name/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(name, 'completed', report['final_radius'], 'energy drift', drift, flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,default=Path('outputs/pockets/balanced-pilot-01'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    args=p.parse_args()
    if args.jobs < 1 or args.jobs > 3:
        p.error('jobs must be 1..3')
    baseline=json.loads((args.baseline/'contract.json').read_text())
    parent=Path(__file__).parent
    for file in ('model.py','pocket_model.h','evolve.cpp'):
        if sha(parent/file) != baseline['sources'][file]:
            raise ValueError(f'Baseline source changed: {file}')
    if sha(EXECUTABLE) != baseline['executable_sha256']:
        raise ValueError('Baseline executable changed')
    modes=radial_modes(20.,2.,2.,.15,.6,31.,h=.01)
    if len(modes) != len(baseline['modes']) or not all(
        m['l']==b['l'] and m['n']==b['n'] and abs(m['omega']-b['omega']) < 1e-12
        for m,b in zip(modes,baseline['modes'])):
        raise ValueError('Baseline modes changed')
    args.output.mkdir(parents=True,exist_ok=False)
    contract=dict(runs=RUNS,baseline=str(args.baseline.resolve()),
                  baseline_contract_sha256=sha(args.baseline/'contract.json'),
                  sources={f:sha(parent/f) for f in ('BALANCED_REFINEMENT.md',
                           'run_balanced_refinement.py','wall_diagnostics.py')},
                  executable_sha256=sha(EXECUTABLE))
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    with ThreadPoolExecutor(args.jobs) as pool:
        reports=dict(zip(RUNS,pool.map(lambda name:execute(name,args.output,baseline,modes),RUNS)))
    base=history(args.baseline/'balanced')
    comparisons={}
    for name in RUNS:
        h=history(args.output/name)
        shared,i,j=np.intersect1d(np.round(base['u'],9),np.round(h['u'],9),return_indices=True)
        dr=float(np.max(np.abs(base['fv_radius'][i]-h['fv_radius'][j])))
        df=float(np.max(np.abs(base['axion_far'][i]-h['axion_far'][j]))/baseline['target_axion_energy'])
        limit=.05 if name=='time' else .25
        comparisons[name]=dict(last_shared_time=float(shared[-1]),max_radius_difference=dr,
                               max_far_storage_difference=df,
                               passed=bool(shared[-1]>=80. and dr<=limit and df<=.005))
    result=dict(runs=reports,comparisons=comparisons,
                numerical_consistency=all(all(r['gates'].values()) for r in reports.values())
                and all(c['passed'] for c in comparisons.values()),
                interpretation='Field-only early consistency study, not continuum or isolated lifetime proof.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(comparisons,indent=2))


if __name__=='__main__':
    main()
