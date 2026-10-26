"""Measured child-process peak memory for identical modest field runs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import numpy as np
from model import radius,config,prepare_run


def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--n',type=int,default=192)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    exes=dict(baseline=args.baseline.resolve(),candidate=args.candidate.resolve())
    (args.output/'contract.json').write_text(json.dumps(dict(n=args.n,executable_sha256={k:digest(v) for k,v in exes.items()},
        source_sha256=digest(Path(__file__)),comparison='Identical fields/parameters, bitwise final snapshot/history required. Child-process RSS only; no large-grid capacity inference.'),indent=2)+'\n')
    r=radius(args.n,16.,[.7,8.1,8.2]);s=np.zeros((args.n,args.n,args.n,4))
    s[...,0]=.5*(1+np.tanh((r-3)/.7));s[...,1]=.15*np.exp(-r*r/6);s[...,3]=.08*np.exp(-r*r/4)
    cfg=config(args.n,16.,.005,.02,lam=8.,dv=.01,mass_tv=2.,mass_power=4,every=2,slices=2,
               monitor_radius=6.,tracking_inner=5.,tracking_outer=7.)
    for label in exes:prepare_run(s,cfg,args.output/label)
    del s,r
    results={}
    for label,exe in exes.items():
        root=args.output/label
        option='-l' if sys.platform=='darwin' else '-v'
        result=subprocess.run(['/usr/bin/time',option,str(exe),'input='+str((root/'input.in').resolve())],capture_output=True,text=True)
        (root/'run.log').write_text(result.stdout+result.stderr);result.check_returncode()
        if sys.platform=='darwin':
            match=re.search(r'(\d+)\s+maximum resident set size',result.stderr);factor=1
        else:
            match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',result.stderr);factor=1024
        if match is None:raise RuntimeError('Cannot parse process peak RSS')
        footprint=re.search(r'(\d+)\s+peak memory footprint',result.stderr)
        results[label]=dict(max_rss_bytes=int(match[1])*factor,
             peak_memory_footprint_bytes=int(footprint[1]) if footprint else None,
             final_sha256=digest(root/'native/final.bin'),history_sha256=digest(root/'native/history.csv'))
        print(label,results[label],flush=True)
    report=dict(runs=results,bitwise_equal=all(results['baseline'][k]==results['candidate'][k] for k in ('final_sha256','history_sha256')),
                interpretation='Modest-grid peak measurement. Read-only mapped file pages may count in RSS; no claim that a larger grid fits.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['bitwise_equal']:raise RuntimeError('Field equivalence failed')


if __name__=='__main__':main()
