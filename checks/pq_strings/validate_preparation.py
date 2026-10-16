"""Independent checks of continuous core-growth to physical-coupling handoff."""
import argparse
import json
from pathlib import Path
import numpy as np
from model import EXECUTABLE,config,run,read_state,reference,radial_coupling,energies

GATES=dict(final_state_error=2e-4,refinement_ratio=3.5,coefficient_error=2e-13,
           snapshot_energy_error=2e-12,relative_energy_balance=.003)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=EXECUTABLE);args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'contract.json').write_text(json.dumps(GATES,indent=2)+'\n')
    checks={};details={}
    n=4;box=4.;coordinates=np.meshgrid((np.arange(n)+.5)*box/n,*([(np.arange(n)+.5)*box/n]*2),indexing='ij')
    s=np.zeros((n,n,n,4));s[...,0]=1.05+.04*np.sin(2*np.pi*coordinates[0]/box)
    s[...,1]=.15+.03*np.cos(2*np.pi*coordinates[1]/box);s[...,2]=.02;s[...,3]=-.03
    for cg in (3.,2.003):
        errors=[]
        for dt in (.005,.0025):
            cfg=config(n,box,dt,2.,lam=.03,origin=2.,star=3.5,exponent=4.,chi_star=.002,every=1)
            cfg['core_growth_until']=cg
            directory=args.output/f'cg{cg}_dt{dt}';cfg=run(s,cfg,directory,args.executable)
            meta,final=read_state(directory/'native/final.bin')
            exact=reference(s,cfg,[2.])[0];error=float(np.max(abs(final-exact)));errors.append(error)
            h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True)
            coupling=np.array([radial_coupling(cfg,u) for u in h['u']])
            coefficient_error=float(np.max(abs(h['lambda']-coupling)))
            e=energies(final,box,meta['a'],radial_coupling(cfg,2.),meta['chi'])
            energy_error=max(abs(e[k]-h[k][-1]) for k in e)
            work=np.trapezoid(h['source']+h['dilution'],h['u'])
            residual=abs(h['total'][-1]-h['total'][0]-work)/max(h['total'][0],h['total'][-1],1e-30)
            checks[f'{cg}_{dt}_state']=error<GATES['final_state_error']
            checks[f'{cg}_{dt}_coefficients']=coefficient_error<GATES['coefficient_error']
            checks[f'{cg}_{dt}_energy']=energy_error<GATES['snapshot_energy_error']
            checks[f'{cg}_{dt}_balance']=residual<GATES['relative_energy_balance']
            checks[f'{cg}_{dt}_source_partition']=np.max(abs(h['source']-h['lambda_source']-h['chi_source']))<2e-14
            details[f'{cg}_{dt}']=dict(final_error=error,coefficient_error=coefficient_error,
                energy_error=float(energy_error),balance=float(residual),
                negative_lambda_work=float(np.trapezoid(h['lambda_source'],h['u'])))
        checks[f'{cg}_refinement']=errors[0]/errors[1]>GATES['refinement_ratio']
        details[f'{cg}_refinement']=errors[0]/errors[1]
    checks={k:bool(v) for k,v in checks.items()}
    report=dict(gates=GATES,checks=checks,details=details,passed=all(checks.values()))
    (args.output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed'] else 1)
if __name__=='__main__':main()
