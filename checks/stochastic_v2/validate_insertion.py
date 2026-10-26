"""Event-aligned coupled fields versus piecewise independent DOP853 and maps."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import read_state
from expanding_model import invoke,expanding_rhs,scale,balances
from stochastic_reference import insertion,energy_totals
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True);p.add_argument('--expansion-validation',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();build=verify_build(args.executable);out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'gates.json').write_text(json.dumps(dict(error=2e-5,ratio=3,jump=1e-10,no_event=1e-12),indent=2)+'\n')
    n=8;box=8.;h0=.2;end=.3;y=np.zeros((n,n,n,4));y[...,1]=.2
    events=np.array([[.05,2,2,2],[.075,2,2,2],[.1,6,6,6.]])
    table=out/'events.txt';np.savetxt(table,events,fmt='%.17g')
    state=y.copy();last=0.;reference_events=[]
    for u,*centre in events:
        sol=solve_ivp(lambda t,v:expanding_rhs(t,v,n,box,h0),(last,u),state.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
        if not sol.success:raise RuntimeError(sol.message)
        before=sol.y[:,-1].reshape(y.shape);state,decision=insertion(before,box,u,h0,centre,1.5,.5)
        reference_events.append(decision['accepted']);last=u
    sol=solve_ivp(lambda t,v:expanding_rhs(t,v,n,box,h0),(last,end),state.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
    if not sol.success:raise RuntimeError(sol.message)
    reference=sol.y[:,-1].reshape(y.shape);np.save(out/'reference.npy',reference);runs=[]
    for dt in (.02,.01,.005):
        dest=out/f'dt-{dt}';sec=invoke(args.executable,dest,y,box,dt,end,h0,cadence=.02,frames=.3,
            extra=dict(candidate_file=str(table),seed_radius=1.5,seed_width=.5,core_tolerance=.05))
        _,final=read_state(dest/'native/final.bin');ledger=np.atleast_1d(np.genfromtxt(dest/'native/births.csv',delimiter=',',names=True))
        h=np.genfromtxt(dest/'native/history.csv',delimiter=',',names=True);jumps=[];maps=[]
        for i,row in enumerate(ledger):
            if row['accepted']:
                m,pre=read_state(dest/f'native/birth_pre_{i+1}.bin');_,post=read_state(dest/f'native/birth_post_{i+1}.bin')
                expected,decision=insertion(pre,box,m['t'],h0,[row['x'],row['y'],row['z']],1.5,.5)
                maps.append(float(np.max(abs(post-expected))))
                jump=energy_totals(post,box,scale(m['t'],h0))-energy_totals(pre,box,scale(m['t'],h0))
                native=np.array([row['delta_total'],row['delta_S'],row['delta_A']])
                jumps.append(float(np.max(abs(jump-native))/max(1,float(np.max(abs(jump))))))
                if not np.array_equal(pre[...,1::2],post[...,1::2]):raise RuntimeError('Insertion changed axion state')
        total,axion=balances(h)
        runs.append(dict(dt=dt,error=float(np.max(abs(final-reference))),accepted=ledger['accepted'].astype(bool).tolist(),
            map_error=max(maps),jump_error=max(jumps),total_balance=float(max(abs(total))/h['total'][0]),
            axion_balance=float(max(abs(axion))/h['total'][0]),seconds=sec))
    ratios=[runs[i]['error']/runs[i+1]['error'] for i in (0,1)]
    # No-event comparison uses exactly the stage-1 initial state and history cadence.
    olddir=args.expansion_validation/'False-0.01';_,initial=read_state(olddir/'native/initial.bin');_,old=read_state(olddir/'native/final.bin')
    dest=out/'no-events';invoke(args.executable,dest,initial,6.,.01,2.,.2,cadence=.01,frames=2.)
    _,new=read_state(dest/'native/final.bin');noevent=float(np.max(abs(new-old)))
    checks=dict(reference=runs[-1]['error']<2e-5,refinement=min(ratios)>3,acceptance=all(r['accepted']==reference_events for r in runs),
        maps=all(r['map_error']<2e-14 for r in runs),jumps=all(r['jump_error']<1e-10 for r in runs),no_event_regression=noevent<1e-12)
    report=dict(passed=bool(all(checks.values())),checks=checks,runs=runs,ratios=ratios,reference_acceptance=reference_events,
        no_event_error=noevent,build=build,script_sha256=sha(__file__),reference_sha256=sha(Path(__file__).with_name('stochastic_reference.py')),
        scope='Small semi-discrete insertion/ODE validation, not resolved supercritical bubble physics')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='build'},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
