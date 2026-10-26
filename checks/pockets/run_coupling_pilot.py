"""Compiled field phi^2/phi^4 comparison; COUPLING_EXPERIMENT.md."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from model import EXECUTABLE, config, eigen_state, energies, history, radial_modes, read_slice, read_state, run
from wall_diagnostics import local_widths


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--validation',type=Path,default=Path('outputs/pockets/validation-coupling-01/validation.json'))
    args=p.parse_args()
    validation=json.loads(args.validation.read_text())
    binary_sha=hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest()
    if not validation['passed'] or validation['executable_sha256']!=binary_sha or not validation['checks']['quartic_ode']:
        raise ValueError('Validated coupling executable required')
    args.output.mkdir(parents=True,exist_ok=False)
    parent=Path(__file__).parent
    contract=dict(runs=dict(p2=dict(mass_power=2),p4=dict(mass_power=4)),
                  N=128,box=64.,r0=20.,width=2.,dt=.0125,duration=80.,target_axion_energy=2680.82573106329,
                  seed=0,omega_band=[.15,.6],radial_h=.01,radial_rmax=31.,r_cut=30.,
                  sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest() for f in (
                      'COUPLING_EXPERIMENT.md','run_coupling_pilot.py','model.py','pocket_model.h','evolve.cpp','wall_diagnostics.py')},
                  executable_sha256=binary_sha,validation_sha256=hashlib.sha256(args.validation.read_bytes()).hexdigest())
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    def execute(power):
        name=f'p{power}';target=contract['target_axion_energy']
        modes=radial_modes(20.,2.,2.,.15,.6,31.,h=.01,mass_power=power)
        s=eigen_state(128,64.,20.,2.,modes,1.,0,30.)
        e=energies(s,64.,2.,.01,4.,mass_power=power)
        amplitude=float(np.sqrt(target/sum(e[k] for k in ('axion_kinetic','axion_gradient','coupling'))))
        s[...,1]*=amplitude;s[...,3]*=amplitude
        state_contract=dict(amplitude=amplitude,mode_count=sum(2*m['l']+1 for m in modes),
                            modes=[dict(l=m['l'],n=m['n'],omega=m['omega']) for m in modes],
                            assumption='Seeded Gaussian classical high-occupation equipartition, not thermalised.')
        (args.output/f'{name}_state.json').write_text(json.dumps(state_contract,indent=2)+'\n')
        cfg=config(128,64.,.0125,80.,lam=2.,dv=.01,mass_tv=2.,every=40,slices=320,mass_power=power)
        error=None
        try:
            run(s,cfg,args.output/name)
        except RuntimeError as exc:
            error=str(exc)
        h=history(args.output/name)
        _,transferred=read_state(args.output/name/'native/initial.bin')
        ax=sum(h[k] for k in ('axion_kinetic','axion_gradient','coupling'))
        drift=float(np.max(np.abs(h['total']/h['total'][0]-1)))
        walls=[]
        for path in sorted((args.output/name/'native').glob('slice_*.bin')):
            meta,plane=read_slice(path)
            walls.append(dict(time=meta['u'],**local_widths(plane[...,0],.5)))
        gates=dict(transfer=bool(np.array_equal(s,transferred)),
                   completed=error is None and bool(np.isclose(h['u'][-1],80.)),energy=drift<5e-4)
        report=dict(gates=gates,solver_error=error,energy_drift=drift,initial_axion=float(ax[0]),
                    radius_loss=float(h['fv_radius'][0]-h['fv_radius'][-1]),final_radius=float(h['fv_radius'][-1]),
                    axion_depletion=float(1-ax[-1]/ax[0]),final_far_storage=float(h['axion_far'][-1]/ax[0]),
                    walls=walls,slice_diagnostic_pass=all(w['diagnostic_pass'] for w in walls),**state_contract)
        (args.output/name/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(name,{k:v for k,v in report.items() if k not in ('walls','modes')},flush=True)
        return report
    with ThreadPoolExecutor(2) as pool:
        outcomes=dict(zip(('p2','p4'),pool.map(execute,(2,4))))
    a,b=outcomes['p2'],outcomes['p4']
    improved=bool(a['radius_loss']>0 and a['axion_depletion']>0 and b['radius_loss']<=.5*a['radius_loss']
                  and b['axion_depletion']<=.5*a['axion_depletion'])
    original=history(Path('outputs/pockets/balanced-pilot-01/balanced'));new=history(args.output/'p2')
    shared,i,j=np.intersect1d(np.round(original['u'],9),np.round(new['u'],9),return_indices=True)
    dr=float(np.max(np.abs(original['fv_radius'][i]-new['fv_radius'][j])))
    df=float(np.max(np.abs(original['axion_far'][i]-new['axion_far'][j]))/contract['target_axion_energy'])
    reproduced=bool(shared[-1]>=80 and dr<.001 and df<1e-5)
    result=dict(runs=outcomes,screening_improvement=improved,baseline_reproduced=reproduced,
                baseline_max_radius_difference=dr,baseline_max_far_difference=df,
                interpretation='Toy coupling-family pilot. Any improvement requires independent refinement and long survival.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Screening improvement:',improved,'Baseline reproduced:',reproduced)


if __name__=='__main__':
    main()
