"""Run predeclared field-only suites with reference-grid amplitude fixing.

Each case retains full raw inputs, build/source identity, and an independent
Python initial-energy calculation. Never rescale the gas on a refined grid.
The JSON contract names a human-readable declaration; no simulation output
directory may exist when starting a new suite.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from streaming_budget import budget
from history_diagnostics import asphericity
from monitor_reference import continuity_summary
import structured_wall
from model import (EXECUTABLE, config, eigen_state, energies, history,
                   lambda_for_width, radial_modes, read_slice, read_state, run)
from model import prepare_run,run_prepared,initial_transfer_exact
from run_prototype import gas_diagnostics
from wall_diagnostics import local_widths


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_case(name, spec, gates, root, executable=EXECUTABLE):
    n,box=spec['n'],spec['box'];power=spec['mass_power']
    lam=lambda_for_width(spec['width']);dv=spec['delta_v'];mass=spec['mass_tv']
    stiffening=spec.get('vacuum_stiffening',0.);profile=None;profile_identity=None
    if stiffening:
        z=np.linspace(-32.,32.,64001);values=structured_wall.kink(lam,stiffening,extent=32.)(z)
        path=root/f'{name}_profile.npz';np.savez(path,z=z,phi=values)
        profile_identity=dict(path=str(path.resolve()),sha256=sha(path),spacing=.001)
        profile=lambda radius:np.interp(np.asarray(radius)-spec['r0'],z,values,left=0.,right=1.)
    modes=radial_modes(spec['r0'],spec['width'],mass,spec['omega_min'],spec['omega_max'],
                       spec['radial_rmax'],h=spec['radial_h'],mass_power=power,profile=profile)
    if not modes:
        raise ValueError(f'No initial modes in {name}')
    unit=eigen_state(spec['normalization_n'],spec['normalization_box'],spec['r0'],spec['width'],
                     modes,1.,spec['seed'],spec['r_cut'],profile=profile)
    energy=energies(unit,spec['normalization_box'],lam,dv,mass**2,mass_power=power,vacuum_stiffening=stiffening)
    amplitude=float(np.sqrt(spec['target_axion_energy']/sum(energy[k] for k in (
        'axion_kinetic','axion_gradient','coupling'))))
    del unit
    s=eigen_state(n,box,spec['r0'],spec['width'],modes,amplitude,spec['seed'],spec['r_cut'],profile=profile)
    independent=energies(s,box,lam,dv,mass**2,mass_power=power,vacuum_stiffening=stiffening)
    initial_axion=sum(independent[k] for k in ('axion_kinetic','axion_gradient','coupling'))
    state_contract=dict(spec=spec,lambda_=lam,amplitude=amplitude,
                        initial_independent_energies=independent,
                        mode_count=sum(2*m['l']+1 for m in modes),
                        modes=[dict(l=m['l'],n=m['n'],omega=m['omega']) for m in modes],
                        profile=profile_identity,
                        occupation='Gaussian classical high-occupation equipartition; no thermalisation',
                        penetration_depth_at_max_frequency=float(1/np.sqrt(mass**2-max(m['omega'] for m in modes)**2)))
    (root/f'{name}_state.json').write_text(json.dumps(state_contract,indent=2,allow_nan=False)+'\n')
    every=round(spec['history_interval']/spec['dt']);slices=round(spec['slice_interval']/spec['dt'])
    if not np.isclose(every*spec['dt'],spec['history_interval']) or not np.isclose(slices*spec['dt'],spec['slice_interval']):
        raise ValueError('Diagnostic intervals must be integral step multiples')
    frames=round(spec.get('snapshot_interval',0.)/spec['dt'])
    cfg=config(n,box,spec['dt'],spec['duration'],lam=lam,dv=dv,mass_tv=mass,
               every=every,slices=slices,frames=frames,mass_power=power,
               monitor_radius=spec.get('monitor_radius',0.),vacuum_stiffening=stiffening,
               tracking_inner=spec.get('tracking_inner',0.),tracking_outer=spec.get('tracking_outer',0.))
    error=None
    prepare_run(s,cfg,root/name)
    del s
    try:
        run_prepared(root/name,executable)
    except RuntimeError as exc:
        error=str(exc)
    return summarize_case(name,spec,gates,root,state_contract,error)


def summarize_case(name,spec,gates,root,state_contract,error=None,output=None):
    n,box=spec['n'],spec['box'];power=spec['mass_power']
    lam=lambda_for_width(spec['width']);dv=spec['delta_v'];mass=spec['mass_tv']
    initial_axion=float(sum(state_contract['initial_independent_energies'][k]
                            for k in ('axion_kinetic','axion_gradient','coupling')))
    amplitude=state_contract['amplitude']
    h=np.atleast_1d(history(root/name))
    transfer=initial_transfer_exact(root/name/'input.bin',root/name/'native/initial.bin')
    ax=sum(h[k] for k in ('axion_kinetic','axion_gradient','coupling'))
    drift=float(np.max(np.abs(h['total']/h['total'][0]-1)))
    walls=[]
    for path in sorted((root/name/'native').glob('slice_*.bin')):
        meta,plane=read_slice(path)
        walls.append(dict(time=meta['u'],**local_widths(plane[...,0],box/n)))
    final=root/name/'native/final.bin';full_wall=None;gas=None;far_dilaton=None
    if final.exists():
        _,state=read_state(final)
        if spec.get('final_3d_wall',False):
            full_wall=local_widths(state[...,0],box/n)
        if np.any(state[...,0]<.5):
            gas=gas_diagnostics(state,box)
        del state
        far_dilaton=budget(final,lam=lam,dv=dv,mass=mass,mass_power=power,vacuum_stiffening=spec.get('vacuum_stiffening',0.))['far_dilaton_energy']
    shape_measure=asphericity(h)
    numerical=dict(transfer=transfer,completed=error is None and bool(np.isclose(h['u'][-1],spec['duration'])),
                   energy=drift<gates['energy_drift'],
                   independent_initial_energy=bool(abs(float(ax[0])-initial_axion)/max(initial_axion,1e-30)<1e-12))
    flow=continuity_summary(h)
    tracking=continuity_summary(h,'tracking')
    if spec.get('monitor_radius',0)>0:
        numerical['monitor_continuity']=flow is not None and max(flow['max_continuity_residual'].values())<=gates['monitor_continuity']
        numerical['monitor_cadence']=flow is not None and max(flow['coarser_cadence_difference'].values())<=gates['monitor_cadence']
    if spec.get('tracking_outer',0)>0:
        numerical['tracking_continuity']=tracking is not None and max(tracking['max_continuity_residual'].values())<=gates['monitor_continuity']
        numerical['tracking_cadence']=tracking is not None and max(tracking['coarser_cadence_difference'].values())<=gates['monitor_cadence']
        numerical['tracking_enclosure']=bool(np.max(h['tracking_max_fv_radius'])<spec['tracking_inner'])
    wall_pass=bool(walls) and all(w['diagnostic_pass'] for w in walls) and (
        full_wall is not None and full_wall['diagnostic_pass'])
    candidate=bool(all(numerical.values()) and wall_pass and h['fv_radius'].min()>=gates['candidate_min_radius']
                   and 1-ax[-1]/ax[0]<=gates['candidate_depletion'])
    report=dict(spec=spec,numerical_gates=numerical,energy_drift=drift,solver_error=error,
                initial_axion_energy=float(ax[0]),final_axion_energy=float(ax[-1]),
                axion_depletion=float(1-ax[-1]/ax[0]),min_radius=float(h['fv_radius'].min()),
                final_radius=float(h['fv_radius'][-1]),last_time=float(h['u'][-1]),
                final_far_axion_storage=float(h['axion_far'][-1]/ax[0]) if h['fv_volume'][-1]>0 else None,
                final_far_dilaton_storage=far_dilaton,
                max_asphericity=float(np.nanmax(shape_measure)) if np.any(np.isfinite(shape_measure)) else None,gas_final=gas,
                walls=walls,full_final_wall=full_wall,wall_diagnostic_pass=bool(wall_pass),
                candidate=candidate,amplitude=amplitude,mode_count=state_contract['mode_count'])
    report['monitor']=flow
    if tracking is not None:
        report['tracking']=tracking
        report['tracking_storage_retention']=float(h['tracking_energy'][-1]/h['tracking_energy'][0])
        report['tracking_minimum_inner_clearance']=float(spec['tracking_inner']-np.max(h['tracking_max_fv_radius']))
    destination=(output or root)/name;destination.mkdir(parents=True,exist_ok=True)
    (destination/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(name,{k:v for k,v in report.items() if k not in ('spec','walls','gas_final','full_final_wall')},flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--suite',type=Path)
    p.add_argument('--analyze-source',type=Path,help='Analyze saved fields into a fresh output; never rerun evolution')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=EXECUTABLE)
    p.add_argument('--jobs',type=int,default=2)
    args=p.parse_args()
    if args.jobs<1 or args.jobs>2:
        p.error('jobs must be 1 or 2')
    if args.analyze_source:
        source=args.analyze_source
        original=json.loads((source/'contract.json').read_text())
        args.output.mkdir(parents=True,exist_ok=False)
        results={}
        for name,spec in original['cases'].items():
            state_contract=json.loads((source/f'{name}_state.json').read_text())
            error=None if (source/name/'native/final.bin').exists() else 'No final snapshot; completion unknown'
            results[name]=summarize_case(name,spec,original['declaration']['gates'],source,state_contract,error,args.output)
        report=dict(cases=results,analysis_source=str(source.resolve()),
                    source_contract_sha256=sha(source/'contract.json'),analysis_sha256=sha(__file__),
                    interpretation='Fresh analysis of existing fields; no physics rerun. Does not establish long survival.')
        (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        return
    if args.suite is None:
        p.error('--suite or --analyze-source required')
    declaration=json.loads(args.suite.read_text());validation_path=Path(declaration['validation'])
    validation=json.loads(validation_path.read_text())
    if not validation['passed'] or validation['executable_sha256']!=sha(args.executable):
        raise ValueError('Current executable must have a passing validation')
    specs={name:dict(declaration['common'],**case) for name,case in declaration['cases'].items()}
    for name in specs:
        if not name or any(ch not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for ch in name):
            raise ValueError('Case names must be simple path components')
    parent=Path(__file__).parent;contract_path=parent/declaration['experiment']
    args.output.mkdir(parents=True,exist_ok=False)
    probe=args.executable.resolve().with_name('pocket_backend_probe')
    runtime=dict(omp_num_threads=os.environ.get('OMP_NUM_THREADS'),omp_proc_bind=os.environ.get('OMP_PROC_BIND'))
    if probe.exists():
        result=subprocess.run([str(probe)],capture_output=True,text=True,check=True)
        (args.output/'backend_probe.log').write_text(result.stdout+result.stderr)
        lines=[line.removeprefix('POCKET_BACKEND ') for line in result.stdout.splitlines() if line.startswith('POCKET_BACKEND ')]
        runtime['probe']=json.loads(lines[-1]);runtime['probe_sha256']=sha(probe)
    build_info=args.executable.resolve().parent/'adapter/build-info.json'
    if build_info.exists():
        runtime['build_info']=json.loads(build_info.read_text())
        runtime['build_info_sha256']=sha(build_info)
    contract=dict(declaration=declaration,cases=specs,suite_sha256=sha(args.suite),
                  sources={f:sha(parent/f) for f in ('model.py','pocket_model.h','pocket_monitor.h','accurate_sum.h','evolve.cpp',
                      'wall_diagnostics.py','field_budget.py','streaming_budget.py','history_diagnostics.py','monitor_reference.py','structured_wall.py','run_prototype.py','run_field_suite.py',
                      'tracking_monitor.h','tracking_monitor_reference.py','pocket_geometry.py','initial_state.h')},
                  experiment_sha256=sha(contract_path),validation_sha256=sha(validation_path),
                  executable=str(args.executable.resolve()),executable_sha256=sha(args.executable),runtime=runtime)
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2,allow_nan=False)+'\n')
    snapshot=args.output/'source_snapshot';snapshot.mkdir()
    for name in contract['sources']:
        (snapshot/name).write_bytes((parent/name).read_bytes())
    with ThreadPoolExecutor(args.jobs) as pool:
        results=dict(zip(specs,pool.map(lambda name:run_case(name,specs[name],declaration['gates'],args.output,args.executable.resolve()),specs)))
    result=dict(cases=results,interpretation='Predeclared field-only suite; candidate does not mean long-lived relic.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
