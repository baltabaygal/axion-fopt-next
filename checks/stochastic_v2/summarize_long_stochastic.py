"""Assess sustained conversion with matched rate and independent refinements."""
import argparse
import json
from pathlib import Path
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('main','space','time-half','output'):p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args();dest=args.output
    if dest.exists():raise FileExistsError(dest)
    paths=dict(main=args.main,space=args.space,time_half=args.time_half)
    reports={k:json.loads((v/'report.json').read_text()) for k,v in paths.items()}
    m=reports['main'];s=reports['space'];t=reports['time_half']
    def start(r):return next((v['start'] for v in r['all_cell_intervals'] if v['duration']>=2-1e-10),None)
    starts={k:start(r) for k,r in reports.items()};space=abs(s['axion_energy_final']/m['axion_energy_final']-1);time=abs(t['axion_energy_final']/m['axion_energy_final']-1)
    shifts={k:abs(starts[k]-starts['main']) if starts[k] is not None and starts['main'] is not None else None for k in ('space','time_half')}
    checks=dict(numerical=all(r['passed'] for r in reports.values()),same_candidates=len({r['candidate_sha256'] for r in reports.values()})==1,
        same_acceptance=all(r['accepted_indices']==m['accepted_indices'] for r in reports.values()),
        completion=all(r['sustained_conversion'] for r in reports.values()),
        spatial_energy=space<.1,temporal_energy=time<.02,
        spatial_onset=shifts['space'] is not None and shifts['space']<=.4,
        temporal_onset=shifts['time_half'] is not None and shifts['time_half']<=.08)
    report=dict(passed=bool(all(checks.values())),checks=checks,starts=starts,onset_shifts=shifts,
        numerical_passed=bool(all(v for k,v in checks.items() if k not in ('completion','spatial_onset','temporal_onset'))),
        goal_reached=bool(all(checks.values())),
        spatial_energy_relative=space,temporal_energy_relative=time,
        source_reports={k:dict(path=str(paths[k]),sha256=sha(paths[k]/'report.json')) for k in paths},
        results=reports,script_sha256=sha(__file__),scope='All-cell S>=.5 sustained>=2 in one realization; not vacuum relaxation or continuum wall accuracy')
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
