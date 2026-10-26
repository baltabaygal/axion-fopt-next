"""Compiled structured-wall checks, plus the entire unchanged default suite."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from model import config,energies,force,history,potential,read_state,reference,run
from monitor_reference import monitor
from validate import small_state


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    gates=dict(force_relative=2e-7,ode_absolute=2e-4,refinement=3.5,
               snapshot_relative=1e-12,continuity_relative=1e-4,hessian_absolute=2e-8)
    parent=Path(__file__).parent
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,
        executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
        sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest() for f in (
            'STRUCTURED_WALL_EXPERIMENT.md','validate_structured.py','model.py','structured_wall.py',
            'pocket_model.h','pocket_monitor.h','accurate_sum.h','evolve.cpp')}),indent=2)+'\n')
    baseline=subprocess.run([sys.executable,str(parent/'validate.py'),'--executable',str(args.executable),
                             '--output',str(args.output/'default')],capture_output=True,text=True)
    (args.output/'default.log').write_text(baseline.stdout+baseline.stderr)
    baseline.check_returncode()
    old=json.loads((args.output/'default/validation.json').read_text())
    checks=dict(default_suite=old['passed']);details={};lam,dv,mass,b=2.,.01,2.,12.
    rng=np.random.default_rng(17);phi=rng.uniform(-.05,1.05,100);a=rng.normal(size=100)*.3;eps=1e-6
    f=force(phi,a,lam,dv,mass**2,4,b)
    fd=((potential(phi+eps,a,lam,dv,mass**2,4,b)-potential(phi-eps,a,lam,dv,mass**2,4,b))/(2*eps),
        (potential(phi,a+eps,lam,dv,mass**2,4,b)-potential(phi,a-eps,lam,dv,mass**2,4,b))/(2*eps))
    error=float(max(np.max(np.abs(x-y)/np.maximum(1.,np.abs(y))) for x,y in zip(f,fd)))
    checks['force']=error<gates['force_relative'];details['force_relative_error']=error
    errors=[];continuities=[]
    s=small_state(6,6.)
    for dt in (.0025,.00125):
        cfg=config(6,6.,dt,1.,lam=lam,dv=dv,mass_tv=mass,mass_power=4,
                   vacuum_stiffening=b,every=1,monitor_radius=2.2)
        dest=args.output/f'small_{dt}';run(s,cfg,dest,args.executable)
        h=history(dest);_,initial=read_state(dest/'native/initial.bin');_,last=read_state(dest/'native/final.bin')
        errors.append(float(np.max(np.abs(last-reference(s,cfg,[0.,1.])[1]))))
        checks[f'{dt}_transfer']=bool(np.array_equal(initial,s))
        energy_errors=[];scale=abs(h['total'][0])
        for row,state in ((0,s),(-1,last)):
            values=energies(state,6.,lam,dv,mass**2,mass_power=4,vacuum_stiffening=b)
            values.update(monitor(state,6.,2.2,lam,dv,mass**2,4,b))
            energy_errors.extend(abs(value-h[key][row])/scale for key,value in values.items())
        checks[f'{dt}_accounting']=max(energy_errors)<gates['snapshot_relative']
        # Finite-difference full Hessian, including mixed curvature.
        fp=force(last[...,0]+eps,last[...,1],lam,dv,mass**2,4,b)
        fm=force(last[...,0]-eps,last[...,1],lam,dv,mass**2,4,b)
        ap=force(last[...,0],last[...,1]+eps,lam,dv,mass**2,4,b)
        am=force(last[...,0],last[...,1]-eps,lam,dv,mass**2,4,b)
        matrix=np.empty(last.shape[:-1]+(2,2));matrix[...,0,0]=(fp[0]-fm[0])/(2*eps)
        matrix[...,0,1]=matrix[...,1,0]=(ap[0]-am[0])/(2*eps);matrix[...,1,1]=(ap[1]-am[1])/(2*eps)
        bound=dt*np.sqrt(12+max(0.,np.linalg.eigvalsh(matrix)[...,-1].max()))
        checks[f'{dt}_hessian']=abs(bound-h['dt_frequency_bound'][-1])<gates['hessian_absolute']
        residuals=[]
        for energy,current in [('monitor_energy',h['monitor_phi_flux']+h['monitor_axion_flux']),
                               ('monitor_phi_energy',h['monitor_phi_flux']+h['monitor_exchange']),
                               ('monitor_axion_energy',h['monitor_axion_flux']-h['monitor_exchange'])]:
            residuals.append(float((h[energy][-1]-h[energy][0]+np.trapezoid(current,h['u']))/scale))
        continuities.append(max(abs(x) for x in residuals))
        details[str(dt)]=dict(ode_error=errors[-1],accounting_error=float(max(energy_errors)),continuity=residuals)
    checks['ode']=max(errors)<gates['ode_absolute'];checks['time_refinement']=errors[0]/errors[1]>gates['refinement']
    checks['continuity']=max(continuities)<gates['continuity_relative']
    checks['continuity_refinement']=continuities[0]/continuities[1]>gates['refinement']
    details['time_refinement']=errors[0]/errors[1];details['continuity_refinement']=continuities[0]/continuities[1]
    for vacuum,curvature in ((0.,13.),(1.,13.12)):
        state=np.zeros((4,4,4,4));state[...,0]=vacuum
        cfg=config(4,4.,.0025,.01,lam=lam,dv=dv,mass_tv=mass,mass_power=4,vacuum_stiffening=b,every=1)
        dest=args.output/f'vacuum_{vacuum}';run(state,cfg,dest,args.executable)
        h=history(dest);_,last=read_state(dest/'native/final.bin')
        checks[f'vacuum_{vacuum}_stationary']=bool(np.array_equal(state,last))
        checks[f'vacuum_{vacuum}_curvature']=abs(h['dt_frequency_bound'][0]-.0025*np.sqrt(12+curvature))<1e-12
    base=(args.output/'small_0.0025/input.in').read_text()
    for label,value in [('negative',-1),('limit',16),('nonfinite','nan')]:
        native=(args.output/f'reject_{label}').resolve()
        text='\n'.join(f'vacuum_stiffening = {value}' if line.startswith('vacuum_stiffening = ') else line for line in base.splitlines())
        text=text.replace(str((args.output/'small_0.0025/native').resolve()),str(native))
        path=args.output/f'bad_{label}.in';path.write_text(text+'\n')
        result=subprocess.run([str(args.executable.resolve()),'input='+str(path.resolve())],capture_output=True,text=True)
        checks[f'reject_{label}']=result.returncode!=0 and not native.exists()
    report=dict(checks={k:bool(v) for k,v in checks.items()},details=details,
                default_check_count=len(old['checks']),passed=bool(all(checks.values())),
                executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest())
    (args.output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed'] else 1)


if __name__=='__main__':
    main()
