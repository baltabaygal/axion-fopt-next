"""Frozen-window comparison; retain failed gates instead of changing thresholds."""
import argparse,json,hashlib,shutil
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--time',type=Path,required=True)
    p.add_argument('--space',type=Path,required=True);p.add_argument('--control',type=Path,required=True)
    p.add_argument('--seeds',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    reports={};sources={}
    for key,root in [('base',a.base),('time',a.time),('space',a.space)]+[(f'seed{i}',s) for i,s in enumerate(a.seeds)]:
        raw=(root/'toy_analysis.json').read_bytes();report=json.loads(raw);reports[key]=report
        sources[key]=dict(path=str(root.resolve()),analysis_sha256=hashlib.sha256(raw).hexdigest(),initial_state_sha256=hashlib.sha256((root/'run/input.bin').read_bytes()).hexdigest(),executable_sha256=json.loads((root/'provenance.json').read_text())['executable_sha256'])
        for name in ('toy_analysis.json','toy_history.csv','contract.json','provenance.json','topology.json'):
            shutil.copyfile(root/name,a.output/f'{key}_{name}')
    baseline=reports['base']['summary']['late_mean_nu_toy']
    relative={key:None if baseline is None or r['summary']['late_mean_nu_toy'] is None else abs(r['summary']['late_mean_nu_toy']/baseline-1) for key,r in reports.items() if key in ('time','space')}
    values=[r['summary']['late_mean_nu_toy'] for k,r in reports.items() if k.startswith('seed')]
    ensemble=None if any(v is None for v in values) else dict(count=len(values),mean=float(np.mean(values)),sample_standard_deviation=float(np.std(values,ddof=1)),min=min(values),max=max(values),description='three finite-box realizations; no cosmological scaling claim')
    control=json.loads((a.control/'topology.json').read_text());shutil.copyfile(a.control/'topology.json',a.output/'control_topology.json')
    checks=dict(resolved_wall_driven_collapse=reports['base']['summary']['gates']['defects_resolved'] and reports['base']['summary']['collapse_bracket'] is not None and any(x['plaquettes']>0 and abs(x['tau']-reports['base']['summary']['collapse_bracket'][1])<1e-9 and x['core_cells']>=4 for x in control),
      temporal_yield_change=relative['time'] is not None and relative['time']<=.01,
      spatial_yield_change=relative['space'] is not None and relative['space']<=.1,
      late_plateau=reports['base']['summary']['gates']['plateau'],
      harmonic_late_regime=reports['base']['summary']['gates']['harmonic_tail'] and reports['base']['summary']['gates']['two_mass_periods'],
      seed_yields_available=ensemble is not None,
      seed_harmonic_windows=all(r['summary']['gates']['two_mass_periods'] and r['summary']['gates']['plateau'] for k,r in reports.items() if k.startswith('seed')))
    for key in ('time','space'):
        reports[key]['summary']['gates']['temporal_refinement']=checks['temporal_yield_change']
        reports[key]['summary']['gates']['spatial_refinement']=checks['spatial_yield_change']
    final=dict(interpretation='Explicitly modified small-grid toy; no paper reproduction',
      qualified_plateau_yield=all(checks.values()),checks={k:bool(v) for k,v in checks.items()},
      base=reports['base']['summary'],time=reports['time']['summary'],space=reports['space']['summary'],
      relative_refinement_changes=relative,seed_ensemble=ensemble,seeds={k:r['summary'] for k,r in reports.items() if k.startswith('seed')},sources=sources)
    (a.output/'comparison.json').write_text(json.dumps(final,indent=2));print(json.dumps(final,indent=2))
if __name__=='__main__':main()
