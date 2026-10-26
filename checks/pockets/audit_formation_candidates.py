"""Zero-lattice-allocation hierarchy and local-backreaction audit."""
import argparse
import hashlib
import json
from pathlib import Path
from formation_hierarchy import candidate,false_branch
from paper_reference import Parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rows=[candidate(dv=dv,n=n) for dv in (.01,.05,.1) for n in (192,256)]
    branches=[false_branch(a,Parameters(dv=.05)) for a in (0.,1.,2.,4.,8.,12.)]
    report=dict(kind=__doc__,candidates=rows,fixed_a_false_vacuum_branch=branches,
        interpretation='Initial hierarchy planning only. No candidate is selected as a stable solution. Thin-wall, isotropic pressure, adiabaticity and long-lived isolation remain unverified.',
        next='Local field scattering and initial-profile preparation, then controlled contracting fields when resources permit.',
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('audit_formation_candidates.py','formation_hierarchy.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(initial_estimate_gates=[dict(dv=r['parameters']['dv'],n=r['n'],gates=r['gates'],rc_over_radius=r['rc_over_radius']) for r in rows],branches=branches),indent=2))

if __name__=='__main__':main()
