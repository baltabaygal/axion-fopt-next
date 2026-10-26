"""Lower bound on visited true-phase cells from audited, replayed full checkpoints."""
import argparse
import json
from pathlib import Path
import numpy as np
from model import read_state
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--checkpoints',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    full=np.genfromtxt(args.source/'main/native/history.csv',delimiter=',',names=True)
    covered=None;rows=[];margin=1e-6
    for folder in args.checkpoints+[args.source]:
        path=folder/'main/native/final.bin';meta,y=read_state(path)
        h=np.genfromtxt(folder/'main/native/history.csv',delimiter=',',names=True)
        prefix=full[full['u']<=meta['t']+1e-10]
        if len(prefix)!=len(h):raise RuntimeError('Checkpoint history schedule differs')
        error=float(max(np.max(abs(h[k]-prefix[k])) for k in h.dtype.names))
        if error>=1e-10:raise RuntimeError('Checkpoint replay mismatch')
        if covered is None:covered=np.zeros(y.shape[:3],dtype=bool)
        if covered.shape!=y.shape[:3]:raise RuntimeError('Checkpoint grid differs')
        below=y[...,0]<.5;previously_true=int((below&covered).sum())
        covered|=y[...,0]>=.5+margin
        rows.append(dict(path=str(path),sha256=sha(path),u=meta['t'],below_now=int(below.sum()),
            unobserved_true=int((~covered).sum()),below_now_previously_observed_true=previously_true,
            replay_error=error,history_sha256=sha(folder/'main/native/history.csv')))
    report=dict(passed=True,all_cells_observed_true=bool(covered.all()),threshold=.5,positive_margin=margin,
        checkpoints=rows,script_sha256=sha(__file__),
        interpretation='N96 replayed numerical trajectory only. A zero unobserved count establishes that every cell exceeded .500001 in at least one saved full checkpoint. It does not imply simultaneous conversion, relaxation, or an exact first-passage time. Refinement controls do not save the same intermediate full checkpoints.')
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
