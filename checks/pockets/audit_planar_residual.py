"""Audit small final true-phase axion energy in saved refinement states.

No solver runs; phase partition is not an asymptotic transmission measurement.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from paper_reference import Parameters,terms


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    suite=json.loads((args.suite/'report.json').read_text());rows=[];p=Parameters(dv=0.)
    for row in suite['rows']:
        root=args.suite/(row['name']+'-packet')/'native';source=root/'final.txt'
        state=np.loadtxt(source);h=np.genfromtxt(root/'history.csv',delimiter=',',names=True)
        phi,a,pp,pa=state.T;dx=512./len(state);va=terms(phi,a,p)[1]
        difference=np.roll(a,-1)-a
        density=.5*pa**2+va+.25*(difference**2+np.roll(difference,1)**2)/dx**2
        rows.append(dict(case=row['name'],n=len(state),dt=row['dt'],amplitude_factor=row['amplitude_factor'],
            true_phase_fraction=float(density[phi>.9].sum()*dx/h['axion_energy'][0]),
            interface_fraction=float(density[(phi>=.1)&(phi<=.9)].sum()*dx/h['axion_energy'][0]),
            history_continuity_relative=row['continuity_relative'],
            state_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    by_name={r['case']:r for r in rows}
    report=dict(kind=__doc__,rows=rows,
        spatial_decrease_ratios=[by_name[a]['true_phase_fraction']/by_name[b]['true_phase_fraction'] for a,b in
            (('time-half','space-mid'),('space-mid','space-fine'))],
        amplitude_normalized_ratio=by_name['space-mid']['true_phase_fraction']/by_name['amplitude-half']['true_phase_fraction'],
        interpretation='These are finite-time phase energies, below the integrated accounting residual. Do not claim a converged tiny transmission rate.',
        suite_report_sha256=hashlib.sha256((args.suite/'report.json').read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
