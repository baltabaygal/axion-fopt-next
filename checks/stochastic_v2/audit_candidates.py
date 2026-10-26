"""Prepare and audit the frozen candidate catalogue, without field evolution."""
import argparse
import json
from pathlib import Path
from stochastic_reference import prepare_candidates
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=prepare_candidates(out)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    sources=[Path(__file__),Path(__file__).with_name('stochastic_reference.py'),Path(__file__).resolve().parents[2]/'src/axion_fopt_next/core.py']
    (out/'provenance.json').write_text(json.dumps({str(s):sha(s) for s in sources},indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
