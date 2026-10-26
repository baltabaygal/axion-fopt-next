"""Prepare small 1D input files and audit them; never launch a simulation."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from planar_initial import PlanarSeed,prepare,audit
from paper_reference import Parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    p=Parameters(dv=0.)
    gates=dict(continuum_left_moving_fraction=1e-25,lattice_left_moving_fraction=1e-6,
        initial_energy_in_tv_fraction=1e-10,initial_energy_near_walls_fraction=1e-10,
        boosted_wall_10_90_cells_minimum=4.)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    rows=[]
    for speed in (0.,.3):
        seed=PlanarSeed(speed=speed)
        duration=220. if speed==0. else 180.
        for n in (2048,4096):
            row=audit(n,seed,p);checks={key:bool(row[key]<value) for key,value in gates.items() if not key.endswith('_minimum')}
            checks['initial_wall_resolution']=row['boosted_wall_10_90_cells']>=gates['boosted_wall_10_90_cells_minimum']
            checks['approximate_image_time']=duration<row['other_wall_signal_to_target_estimate']
            paths=[]
            for control in ('packet','empty-wall'):
                directory=args.output/f'u{speed}-n{n}-{control}';directory.mkdir()
                x,state=prepare(n,seed,p)
                if control=='empty-wall':state[:,1]=0.;state[:,3]=0.
                np.savetxt(directory/'input.txt',state,fmt='%.17g')
                params=dict(N=n,lSide=seed.box,dt=.025,tMax=duration,expansion=False,evolver='PV2',PS_type=2,
                    **{'lambda':p.lam},delta_v=p.dv,mass_fv=p.mass_false,mass_kappa=p.kappa,axion_fa=p.fa,
                    output_every=20,initial_state=str((directory/'input.txt').resolve()),
                    result_dir=str((directory/'native').resolve()))
                (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
                (directory/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
                paths.append(str(directory))
            row.update(checks=checks,passed=all(checks.values()),prepared_cases=paths,
                prospective_duration=duration,prospective_dt=.025,additional_checks='Time refinement, larger box, half amplitude and local flux diagnostics still required; do not treat this pair as full validation')
            rows.append(row)
    names=['prepare_planar_scattering.py','planar_initial.py','paper_reference.py']
    report=dict(kind=__doc__,solver_launched=False,rows=rows,gates=gates,passed=all(r['passed'] for r in rows),
        classical_assumption='Deterministic classical weak-field test. Quantum occupations and a physical field-unit map have not been inferred.',
        spectrum='Positive-k deterministic Gaussian coefficients, fixed continuum frequencies for momentum, 6-sigma truncation; no random phases or seed',
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],rows=[{k:r[k] for k in ('n','seed','modes','omega_band','continuum_left_moving_fraction',
        'lattice_left_moving_fraction','initial_energy_in_tv_fraction','initial_energy_near_walls_fraction',
        'initial_modal_energy_above_wall_threshold','continuum_kink_discrete_residual','centre_arrival_estimate',
        'other_wall_signal_to_target_estimate','passed')} for r in rows]),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
