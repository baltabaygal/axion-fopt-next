"""Exact outputs across diagnostic-memory implementations; fields unchanged."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from model import radius,config,run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--n',type=int,default=48)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executables=dict(baseline=args.baseline.resolve(),candidate=args.candidate.resolve())
    parent=Path(__file__).parent
    (args.output/'contract.json').write_text(json.dumps(dict(executables={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in executables.items()},n=args.n,
        sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest() for f in ('verify_memory_backend.py','model.py','evolve.cpp','pocket_model.h','pocket_monitor.h','tracking_monitor.h','initial_state.h')},
        gate='Bitwise full/slice binary and history equality; includes periodic-cut moving pocket, p=2/4,b=0/12.'),indent=2)+'\n')
    r=radius(args.n,16.,[.7,8.1,8.2]);s=np.zeros((args.n,args.n,args.n,4))
    s[...,0]=.5*(1+np.tanh((r-3.)/.7));s[...,1]=.15*np.exp(-r*r/6);s[...,3]=.08*np.exp(-r*r/4)
    s[...,2]=-.12*(np.roll(s[...,0],-1,0)-np.roll(s[...,0],1,0))/(2*(16/args.n))
    s[...,2]+=.001*np.random.default_rng(88).normal(size=r.shape)
    checks={};timings={}
    for power in (2,4):
      for b in (0.,12.):
        name=f'p{power}_b{int(b)}';paths=[];timings[name]={}
        cfg=config(args.n,16.,.005,.2,lam=8.,dv=.01,mass_tv=2.,mass_power=power,vacuum_stiffening=b,
                   every=4,frames=20,slices=10,monitor_radius=6.,tracking_inner=5.,tracking_outer=7.)
        for label,exe in executables.items():
            path=args.output/f'{name}_{label}';start=time.perf_counter();run(s,cfg,path,exe)
            timings[name][label]=time.perf_counter()-start;paths.append(path/'native')
        names=sorted(p.name for p in paths[0].glob('*.bin'))+['history.csv']
        checks[name]={filename:paths[0].joinpath(filename).read_bytes()==paths[1].joinpath(filename).read_bytes() for filename in names}
        checks[name]['same_files']=sorted(p.name for p in paths[0].glob('*.bin'))==sorted(p.name for p in paths[1].glob('*.bin'))
        print(name,all(checks[name].values()),timings[name],flush=True)
    report=dict(checks=checks,timings_seconds=timings,passed=all(all(v.values()) for v in checks.values()),
                interpretation='Exact numerical equivalence of memory paths. Modest wall timings include I/O/startup, not large-grid scaling or peak memory.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:raise RuntimeError('Exact memory-backend comparison failed')


if __name__=='__main__':main()
