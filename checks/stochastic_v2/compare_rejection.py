"""Compare the accelerated extension with its ordinary saved u20 control."""
import argparse
import json
from pathlib import Path
import numpy as np
from model import read_state
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ordinary',type=Path,required=True);p.add_argument('--accelerated',type=Path,required=True)
    p.add_argument('--output',type=Path);args=p.parse_args()
    target=args.output or args.accelerated/'acceleration-audit.json'
    if target.exists():raise FileExistsError(target)
    a=json.loads((args.ordinary/'report.json').read_text());b=json.loads((args.accelerated/'report.json').read_text())
    _,x=read_state(args.ordinary/'main/native/final.bin');_,y=read_state(args.accelerated/'main/native/final.bin')
    error=float(np.max(abs(x-y)));energy=abs(b['axion_energy_final']/a['axion_energy_final']-1)
    checks=dict(local=a['passed'] and b['passed'],state=error<.002,energy=energy<.001,
        accepted=a['accepted_indices']==b['accepted_indices'],prefix=b['prefix_overlap_error']==0)
    report=dict(passed=bool(all(checks.values())),checks=checks,state_error=error,component_errors=np.max(abs(x-y),axis=(0,1,2)).tolist(),
        energy_relative=energy,ordinary_sha256=sha(args.ordinary/'report.json'),accelerated_sha256=sha(args.accelerated/'report.json'),script_sha256=sha(__file__))
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
