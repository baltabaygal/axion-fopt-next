"""Bounded planar survival screen with passive two-wall energy accounting."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from planar_distributed_initial import prepare,SlabSeed
from validate_slab_monitor import execute,integral,reference
from paper_reference import Parameters

HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--executable',type=Path,required=True)
    ap.add_argument('--n',type=int,default=2048);ap.add_argument('--dt',type=float,default=.0125)
    ap.add_argument('--duration',type=float,default=1280.);ap.add_argument('--cadence',type=float,default=.25)
    ap.add_argument('--analyze-existing',action='store_true')
    ap.add_argument('--box',type=float,default=512.)
    ap.add_argument('--monitor-validation',type=Path,default=HERE.parent.parent/'outputs/pockets/slab-monitor-validation-02/report.json')
    a=ap.parse_args()
    if not 128<=a.n<=16384 or any(not np.isfinite(v) or v<=0 for v in (a.dt,a.duration,a.cadence,a.box)) or any(abs(v/a.dt-round(v/a.dt))>1e-8 for v in (a.duration,a.cadence)):raise ValueError('Invalid grid/time contract')
    exe=a.executable.resolve();build=json.loads((exe.parent/'build.json').read_text())
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=build['executable_sha256']:raise ValueError('Binary hash mismatch')
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Stale native build')
    validation=a.monitor_validation
    valid=json.loads(validation.read_text())
    if not valid['passed'] or valid['build']['executable_sha256']!=build['executable_sha256']:raise ValueError('Require matching monitor validation')
    if not a.analyze_existing:a.output.mkdir(parents=True,exist_ok=False)
    elif not (a.output/'contract.json').exists() or (a.output/'report.json').exists():raise ValueError('Require unfinished saved screen')
    gates=dict(energy_over_initial_axion=.002,exchange_over_initial_axion=.002,moving_balance_over_initial_axion=.002,width_cells=4.,endpoint_reference=1e-11)
    contract=dict(n=a.n,box=a.box,dt=a.dt,duration=a.duration,cadence=a.cadence,padding=4.,spatial_order=4,gates=gates,
        time='Minkowski, a=1',units='dimensionless energy per transverse area',solver_timeout_seconds=30,
        hypothesis='Distributed trapped fields may arrest contraction and survive repeated interactions',
        interpretation='Finite planar screen only; neither 3D stability nor irreversible escape',
        validation_report_sha256=hashlib.sha256(validation.read_bytes()).hexdigest())
    if not a.analyze_existing:
        (a.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
        (a.output/'executed_runner.py').write_bytes(Path(__file__).read_bytes())
    elif json.loads((a.output/'contract.json').read_text())!=contract:raise ValueError('Saved contract mismatch')
    _,s,modes=prepare(a.n,SlabSeed(box=a.box));started=time.monotonic()
    if not a.analyze_existing:
        h,f=execute(exe,a.output/'filled',s,a.box,a.dt,a.duration,4.,4,a.cadence)
        elapsed=time.monotonic()-started
        (a.output/'timing.json').write_text(json.dumps(dict(solver_and_io_seconds=elapsed)))
    else:
        h=np.genfromtxt(a.output/'filled/native/history.csv',delimiter=',',names=True)
        f=np.loadtxt(a.output/'filled/native/final.txt')
        elapsed=json.loads((a.output/'timing.json').read_text())['solver_and_io_seconds'] if (a.output/'timing.json').exists() else None
    t=h['t'];e0=h['axion_energy'][0]
    energy=float(np.max(abs(h['total']-h['total'][0]))/e0)
    exchange=float(np.max(abs(h['axion_energy']-e0-integral(h['exchange'],t)))/e0)
    moving=float(np.max(abs(h['slab_axion_energy']-h['slab_axion_energy'][0]+integral(h['slab_axion_flux'],t)-integral(h['slab_exchange'],t)))/e0)
    agreement=max(abs(v-h['slab_'+k][i]) for i,state in ((0,s),(-1,f)) for k,v in reference(state,a.box,Parameters(dv=.0002),4.,4).items())
    length=h['slab_right']-h['slab_left'];rate=h['slab_right_velocity']-h['slab_left_velocity']
    turns=np.flatnonzero((rate[:-1]<0)&(rate[1:]>=0))+1
    outside=(h['axion_energy']-h['slab_axion_energy'])/h['axion_energy']
    checks=dict(energy=energy<.002,exchange=exchange<.002,moving_balance=moving<.002,
        resolved=bool(h['slab_min_width_cells'].min()>=4),endpoint_reference=bool(agreement<1e-11))
    result=dict(analysis_recovered=a.analyze_existing,executed_runner_sha256=hashlib.sha256((a.output/'executed_runner.py').read_bytes()).hexdigest(),contract=contract,build=build,modes=modes,checks=checks,numerical_gates_passed=all(checks.values()),
        solver_and_io_seconds=elapsed,energy_drift_over_initial_axion=energy,exchange_residual=exchange,moving_balance_residual=moving,
        endpoint_reference_absolute=float(agreement),minimum_length=float(length.min()),final_length=float(length[-1]),
        final_left_velocity=float(h['slab_left_velocity'][-1]),final_right_velocity=float(h['slab_right_velocity'][-1]),
        sampled_contraction_to_expansion_times=t[turns].tolist(),final_outside_energy_fraction=float(outside[-1]),
        maximum_outside_energy_fraction=float(outside.max()),final_axion_energy_ratio=float(h['axion_energy'][-1]/e0),
        integrated_net_outward_moving_flux_over_initial=float(integral(h['slab_axion_flux'],t)[-1]/e0),
        integrated_lab_flux_over_initial=float(integral(h['slab_axion_lab_flux'],t)[-1]/e0),
        integrated_advection_over_initial=float(integral(h['slab_axion_advection'],t)[-1]/e0),
        minimum_width_cells=float(h['slab_min_width_cells'].min()),survival_established=False,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_slab_survival.py','validate_slab_monitor.py','slab_retention.py','planar_distributed_initial.py','paper_reference.py','planar_stencils.py','planar_diagnostics.py')})
    (a.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('build','modes','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
