"""Independent expanding coupled-field validation and exact static regression."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import read_state
from expanding_model import expanding_rhs,invoke,physical_terms,balances,scale,cosmic
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--static-reference',type=Path,required=True);args=p.parse_args();build=verify_build(args.executable)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'gates.json').write_text(json.dumps(dict(error=2e-5,ratio=3.5,energy=1e-11,balance=1e-4),indent=2)+'\n')
    cases=[]
    for homogeneous in (True,False):
        n=4;box=6.;h0=.2;y=np.zeros((n,n,n,4));y[:]=[.4,.17,.03,-.02]
        x=(np.arange(n)+.5)/n*2*np.pi
        if not homogeneous:y[...,0]+=.07*np.cos(x)[:,None,None];y[...,1]+=.03*np.sin(x)[None,:,None];y[...,2]+=.02*np.cos(x)[None,None,:]
        ref=solve_ivp(lambda t,v:expanding_rhs(t,v,n,box,h0),(0,2),y.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
        if not ref.success:raise RuntimeError(ref.message)
        np.save(out/f'reference-{homogeneous}.npy',ref.y[:,-1].reshape(y.shape));runs=[]
        for dt in (.02,.01,.005):
            dest=out/f'{homogeneous}-{dt}';seconds=invoke(args.executable,dest,y,box,dt,2,h0,cadence=dt,frames=2)
            _,initial=read_state(dest/'native/initial.bin');_,final=read_state(dest/'native/final.bin');h=np.genfromtxt(dest/'native/history.csv',delimiter=',',names=True)
            energy=sum(float(v.sum()) for v in physical_terms(final,box/n,1.4).values())*(box/n)**3
            total,axion=balances(h)
            runs.append(dict(dt=dt,error=float(np.max(abs(final-ref.y[:,-1].reshape(y.shape)))),
                initial_exact=bool(np.array_equal(initial,y)),snapshot_energy=abs(energy-h['total'][-1])/h['total'][0],
                total_balance=float(np.max(abs(total))/h['total'][0]),axion_balance=float(np.max(abs(axion))/h['total'][0]),
                scale_error=float(np.max(abs(h['a']-scale(h['u'],h0)))),cosmic_error=float(np.max(abs(h['t_cosmic']-cosmic(h['u'],h0)))),seconds=seconds))
        ratios=[runs[i]['error']/runs[i+1]['error'] for i in (0,1)]
        passed=runs[-1]['error']<2e-5 and min(ratios)>3.5 and max(runs[-1]['total_balance'],runs[-1]['axion_balance'])<1e-4
        passed=passed and all(r['initial_exact'] and r['snapshot_energy']<1e-11 and r['scale_error']<1e-12 and r['cosmic_error']<1e-12 for r in runs)
        cases.append(dict(homogeneous=homogeneous,runs=runs,ratios=ratios,passed=bool(passed)))
    regression=[]
    for folder in args.static_reference.iterdir():
        if not folder.is_dir() or not (folder/'native/initial.bin').exists():continue
        _,initial=read_state(folder/'native/initial.bin');_,old=read_state(folder/'native/final.bin');dt=float(folder.name.split('-')[-1])
        dest=out/('static-'+folder.name);invoke(args.executable,dest,initial,6.,dt,.4,0.,cadence=.02,frames=.4)
        _,new=read_state(dest/'native/final.bin');regression.append(dict(case=folder.name,exact=bool(np.array_equal(old,new))))
    report=dict(passed=bool(all(c['passed'] for c in cases) and regression and all(r['exact'] for r in regression)),cases=cases,regression=regression,build=build,
        script_sha256=sha(__file__),reference_sha256=sha(Path(__file__).with_name('expanding_model.py')))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='build'},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
