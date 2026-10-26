"""Append independent Poisson time intervals without changing the saved past."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from stochastic_reference import rate,sample_candidates
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,required=True);p.add_argument('--end',type=float,choices=[20.,22.,24.,26.,28.],required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    original=np.atleast_2d(np.loadtxt(args.original));parts=[original];audits=[]
    for i,(start,end) in enumerate(zip([12.,16.,18.,20.,22.,24.,26.],[16.,18.,20.,22.,24.,26.,28.])):
        if end>args.end:break
        u=np.linspace(start,end,int(np.ceil((end-start)/(12/1024)))+1);q=rate(u)
        seed=np.random.SeedSequence([0,20261020,i])
        c=sample_candidates(u,q,32.,1.,seed,max_candidates=1000000)
        parts.append(np.column_stack((c.times,c.positions)))
        np.savetxt(out/f'rate-{i}.csv',np.column_stack((u,q)),delimiter=',',header='u,q_Gamma_a4',comments='')
        integral=float(np.trapezoid(q,u));truth=quad(rate,start,end,epsabs=1e-12,epsrel=1e-12)[0]
        error=abs(integral/truth-1);assert error<.001
        audits.append(dict(start=start,end=end,seed_entropy=[0,20261020,i],nodes=len(u),expected_count=32**3*integral,
            count=c.candidate_count,quadrature_error=error))
    data=np.concatenate(parts)
    # Retain the original bytes as the literal prefix, including its header.
    target=out/'candidates.txt';target.write_bytes(args.original.read_bytes())
    with target.open('a') as stream:np.savetxt(stream,np.concatenate(parts[1:]),fmt='%.17g')
    reread=np.loadtxt(target);assert np.array_equal(reread[:len(original)],original)
    report=dict(passed=True,end=args.end,box=32.,H0=.025,beta=.4,prefix_count=len(original),candidate_count=len(data),
        original_sha256=sha(args.original),candidate_sha256=sha(target),intervals=audits,numpy=np.__version__,
        script_sha256=sha(__file__),reference_sha256=sha(Path(__file__).with_name('stochastic_reference.py')))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
