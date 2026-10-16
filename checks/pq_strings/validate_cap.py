"""Independent semi-discrete reference for the explicitly modified toy mass law."""
import argparse,json
from pathlib import Path
import numpy as np
from model import config,run,read_state,reference,coefficients

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--executable',type=Path);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    gates=dict(max_error=2e-5,min_refinement_ratio=3.5,coefficient_error=2e-14,source_balance=.001)
    (a.output/'contract.json').write_text(json.dumps(gates,indent=2))
    s=np.zeros((4,4,4,4));s[...,0]=1.03;s[...,1]=.2;s[...,3]=.02
    s[...,0]+=.03*np.cos(2*np.pi*(np.arange(4)+.5)/4)[:,None,None]
    errors=[];checks={};details=[]
    for dt in (.01,.005):
        cfg=config(4,4,dt,4,origin=2,star=3,lam=.08,every=1)
        cfg.update(mass_schedule='smooth_cap',mass_onset=2.703,mass_ramp_duration=1.31,mass_cap=.15)
        root=a.output/f'dt{dt}';cfg=run(s,cfg,root,**({} if a.executable is None else dict(executable=a.executable)))
        _,final=read_state(root/'native/final.bin');err=float(np.max(abs(final-reference(s,cfg,[4])[0])));errors.append(err)
        h=np.genfromtxt(root/'native/history.csv',delimiter=',',names=True)
        ch=np.array([coefficients(cfg,u)[2] for u in h['u']])
        residual=abs(h['total'][-1]-h['total'][0]-np.trapezoid(h['source']+h['dilution'],h['u']))/h['total'][0]
        checks[f'{dt}_state']=err<gates['max_error'];checks[f'{dt}_chi']=np.max(abs(ch-h['chi']))<gates['coefficient_error']
        checks[f'{dt}_balance']=residual<gates['source_balance']
        details.append(dict(dt=dt,error=err,energy_balance=float(residual)))
    checks['refinement']=errors[0]/errors[1]>gates['min_refinement_ratio']
    report=dict(checks={k:bool(v) for k,v in checks.items()},details=details,ratio=errors[0]/errors[1],passed=all(checks.values()))
    (a.output/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
    raise SystemExit(0 if report['passed'] else 1)
if __name__=='__main__':main()
