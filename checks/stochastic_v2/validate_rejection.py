"""Exercise certified rejections versus an independent homogeneous ODE."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import read_state
from expanding_model import invoke,expanding_rhs
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--insertion-validation',type=Path,required=True);args=p.parse_args()
    build=verify_build(args.executable);validation=json.loads(args.insertion_validation.read_text())
    assert validation['passed'] and validation['build']['executable_sha256']==build['executable_sha256']
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    y=np.zeros((4,4,4,4));y[:]=[1,.17,.03,-.02]
    reference=solve_ivp(lambda t,v:expanding_rhs(t,v,4,8.,.2),(0,2.),y.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
    assert reference.success
    events=np.column_stack((np.linspace(.013,1.993,151),np.tile([3.,3.,3.],(151,1))));np.savetxt(out/'events.txt',events,fmt='%.17g')
    runs=[]
    for dt in [.02,.01,.005]:
        dest=out/f'dt-{dt}'
        invoke(args.executable,dest,y,8.,dt,2.,.2,cadence=.02,frames=2.,extra=dict(candidate_file=str(out/'events.txt'),seed_radius=1.5,seed_width=.5,core_tolerance=.05,skip_rejected_after=0))
        _,final=read_state(dest/'native/final.bin');ledger=np.genfromtxt(dest/'native/births.csv',delimiter=',',names=True)
        runs.append(dict(dt=dt,error=float(np.max(abs(final-reference.y[:,-1].reshape(y.shape)))),
            certified=int(ledger['certified_rejection'].sum()),accepted=int(ledger['accepted'].sum()),
            minimum_certificate=float(ledger['witness_lower_bound'][ledger['certified_rejection']>0].min())))
    ratios=[runs[i]['error']/runs[i+1]['error'] for i in [0,1]]
    checks=dict(reference=runs[-1]['error']<2e-5,refinement=min(ratios)>3.5,rejected=all(r['accepted']==0 for r in runs),
        certified=all(r['certified']>0 and r['minimum_certificate']>.05 for r in runs))
    global_runs=[]
    for dt in [.02,.01,.005]:
        dest=out/f'global-{dt}'
        invoke(args.executable,dest,y,4.,dt,2.,.2,cadence=.02,frames=2.,extra=dict(seed_radius=1.5,seed_width=.1,core_tolerance=.05,skip_rejected_after=0,certify_empty_after=0))
        _,final=read_state(dest/'native/final.bin');proof=np.genfromtxt(dest/'native/empty_eligibility.csv',delimiter=',',names=True)
        global_runs.append(dict(dt=dt,error=float(np.max(abs(final-reference.y[:,-1].reshape(y.shape)))),
            steps=len(proof),minimum_witness=float(proof['min_witness'].min()),coverage=bool(abs(proof['u_start'][0])<1e-12 and
            abs(proof['u_end'][-1]-2)<1e-10 and np.max(abs(proof['u_start'][1:]-proof['u_end'][:-1]))<1e-12)))
    global_ratios=[global_runs[i]['error']/global_runs[i+1]['error'] for i in [0,1]]
    checks['global_reference']=global_runs[-1]['error']<2e-5 and min(global_ratios)>3.5
    checks['global_coverage']=all(r['coverage'] and r['minimum_witness']>.05 and r['steps']>0 for r in global_runs)
    dest=out/'stop-control'
    invoke(args.executable,dest,y,4.,.02,2.,.2,cadence=.04,frames=2.,extra=dict(seed_radius=1.5,seed_width=.1,skip_rejected_after=0,certify_empty_after=0,stop_sustained=.12))
    meta,stopped=read_state(dest/'native/final.bin')
    stopref=solve_ivp(lambda t,v:expanding_rhs(t,v,4,4.,.2),(0,.12),y.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
    stoperror=float(np.max(abs(stopped-stopref.y[:,-1].reshape(y.shape))))
    checks['stop_time']=abs(meta['t']-.12)<1e-12 and stoperror<2e-5
    report=dict(passed=bool(all(checks.values())),checks=checks,runs=runs,ratios=ratios,global_runs=global_runs,global_ratios=global_ratios,
        stop_control=dict(time=meta['t'],error=stoperror),build=build,
        insertion_validation_sha256=sha(args.insertion_validation),script_sha256=sha(__file__))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='build'},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
