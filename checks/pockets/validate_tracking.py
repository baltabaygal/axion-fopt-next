"""Compiled moving-window current against its independent reference."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from model import config,history,laplacian,radius,read_state,reference,run
from tracking_monitor_reference import tracking_monitor


def state():
    n=12;box=8.;dx=box/n;r=radius(n,box,[3.4,4.2,4.6]);s=np.zeros((n,n,n,4))
    s[...,0]=.5*(1+np.tanh((r-1.6)/.5))
    s[...,1]=.12*np.exp(-r*r/1.5);s[...,3]=.08*np.exp(-r*r/2)
    s[...,2]=-.1*(np.roll(s[...,0],-1,0)-np.roll(s[...,0],1,0))/(2*dx)
    s[...,2]+=.01*np.random.default_rng(23).normal(size=(n,n,n))
    return s


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    parent=Path(__file__).parent
    gates=dict(accounting_relative=1e-12,geometry_absolute=1e-12,ode_absolute=2e-4,
               refinement=3.5,continuity_relative=1e-4)
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,
        executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
        sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest() for f in (
            'TRACKING_MONITOR_CONTRACT.md','validate_tracking.py','tracking_monitor.h',
            'tracking_monitor_reference.py','pocket_geometry.py','model.py','evolve.cpp','initial_state.h')}),indent=2)+'\n')
    baseline=subprocess.run([sys.executable,str(parent/'validate_structured.py'),'--executable',str(args.executable),
                             '--output',str(args.output/'structured')],capture_output=True,text=True)
    (args.output/'structured.log').write_text(baseline.stdout+baseline.stderr);baseline.check_returncode()
    old=json.loads((args.output/'structured/validation.json').read_text())
    checks=dict(structured_and_default=old['passed']);details={}
    s=state()
    for stiffening in (0.,12.):
        errors=[];continuities=[]
        for dt in (.005,.0025):
            cfg=config(12,8.,dt,.2,lam=2.,dv=.01,mass_tv=2.,mass_power=4,
                       vacuum_stiffening=stiffening,every=1,tracking_inner=2.5,tracking_outer=3.5)
            dest=args.output/f'b{stiffening}_{dt}';run(s,cfg,dest,args.executable)
            h=history(dest);_,last=read_state(dest/'native/final.bin')
            _,initial=read_state(dest/'native/initial.bin')
            checks[f'{stiffening}_{dt}_transfer']=bool(np.array_equal(initial,s))
            errors.append(float(np.max(np.abs(last-reference(s,cfg,[0.,.2])[1]))))
            scale=abs(h['tracking_energy'][0]);accounting=[];geometry=[]
            for row,value in ((0,s),(-1,last)):
                expected=tracking_monitor(value,8.,2.5,3.5,2.,.01,4.,4,stiffening)
                for key in ('tracking_energy','tracking_phi_energy','tracking_axion_energy','tracking_phi_flux',
                            'tracking_axion_flux','tracking_exchange','tracking_lab_phi_flux','tracking_lab_axion_flux',
                            'tracking_phi_advection','tracking_axion_advection'):
                    accounting.append(float(abs(expected[key]-h[key][row])/scale))
                for j,axis in enumerate('xyz'):
                    geometry.append(abs(expected['centre'][j]-h['tracking_c'+axis][row]))
                    geometry.append(abs(expected['velocity'][j]-h['tracking_v'+axis][row]))
            checks[f'{stiffening}_{dt}_accounting']=max(accounting)<gates['accounting_relative']
            checks[f'{stiffening}_{dt}_geometry']=max(geometry)<gates['geometry_absolute']
            residuals=[]
            for energy,current in [('tracking_energy',h['tracking_phi_flux']+h['tracking_axion_flux']),
                ('tracking_phi_energy',h['tracking_phi_flux']+h['tracking_exchange']),
                ('tracking_axion_energy',h['tracking_axion_flux']-h['tracking_exchange'])]:
                residuals.append(float((h[energy][-1]-h[energy][0]+np.trapezoid(current,h['u']))/scale))
            continuities.append(max(abs(x) for x in residuals))
            details[f'{stiffening}_{dt}']=dict(ode_error=errors[-1],accounting=max(accounting),geometry=float(max(geometry)),continuity=residuals)
        checks[f'{stiffening}_ode']=max(errors)<gates['ode_absolute']
        checks[f'{stiffening}_time_refinement']=errors[0]/errors[1]>gates['refinement']
        checks[f'{stiffening}_continuity']=max(continuities)<gates['continuity_relative']
        checks[f'{stiffening}_continuity_refinement']=continuities[0]/continuities[1]>gates['refinement']
        details[f'{stiffening}_refinement']=dict(state=errors[0]/errors[1],continuity=continuities[0]/continuities[1])
    base=(args.output/'b0.0_0.005/input.in').read_text()
    for label,key,value in [('negative','tracking_inner',-1),('order','tracking_inner',3.6),
                             ('large','tracking_outer',4),('nonfinite','tracking_inner','nan')]:
        native=(args.output/f'reject_{label}').resolve()
        text='\n'.join(f'{key} = {value}' if line.startswith(key+' = ') else line for line in base.splitlines())
        text=text.replace(str((args.output/'b0.0_0.005/native').resolve()),str(native))
        path=args.output/f'bad_{label}.in';path.write_text(text+'\n')
        result=subprocess.run([str(args.executable.resolve()),'input='+str(path.resolve())],capture_output=True,text=True)
        checks[f'reject_{label}']=result.returncode!=0 and not native.exists()
    report=dict(checks={k:bool(v) for k,v in checks.items()},details=details,passed=bool(all(checks.values())),
                executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                interpretation='Static moving-window numerical validation, not physical pocket survival.')
    (args.output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed'] else 1)


if __name__=='__main__':
    main()
