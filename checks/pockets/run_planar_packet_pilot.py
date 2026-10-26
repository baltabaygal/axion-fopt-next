"""One bounded stationary-wall packet/control pair; preliminary, not converged."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from paper_reference import Parameters,terms
from planar_diagnostics import wall_geometry

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--prepared',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--detector-design',type=Path,help='Passing continuum detector-design report; default retains original pilot geometry')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve();build=json.loads((executable.parent/'build.json').read_text())
    for name,digest in build['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale line build')
    if hashlib.sha256(executable.read_bytes()).hexdigest()!=build['executable_sha256']:raise RuntimeError('Wrong binary')
    gates=dict(maximum_total_solver_seconds=45.,axion_continuity_relative=2e-3,
        output_cadence_relative=2e-3,incident_crossing_fraction_min=.99,final_window_missing_fraction_max=.02,
        empty_axion_exact=True)
    window=dict(inner=80.,outer=112.,centre=96.,velocity=0.)
    design_identity=None
    if args.detector_design is not None:
        design=json.loads(args.detector_design.read_text())
        if not design['passed']:raise ValueError('Detector reference design did not pass')
        if design['contract']['n']!=2048 or design['contract']['box']!=512. or design['contract']['time_end']!=220.:
            raise ValueError('Detector design does not match pilot physical configuration')
        window.update(design['contract']['proposed'])
        design_identity=dict(path=str(args.detector_design.resolve()),sha256=hashlib.sha256(args.detector_design.read_bytes()).hexdigest())
    contract=dict(gates=gates,n=2048,box=512.,duration=220.,dt=.025,output_interval=.5,
        window=window,detector_design=design_identity,
        scope='Single-grid stationary weak-field pilot; no spatial, temporal, amplitude or boundary-convergence claim')
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    environment=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    results={};spent=0.
    for label in ('empty-wall','packet'):
        source=args.prepared/f'u0.0-n2048-{label}';params=json.loads((source/'parameters.json').read_text())
        directory=args.output/label;directory.mkdir();initial=np.loadtxt(source/'input.txt');np.savetxt(directory/'input.txt',initial,fmt='%.17g')
        params.update(line_window_inner=window['inner'],line_window_outer=window['outer'],line_window_centre=window['centre'],line_window_velocity=window['velocity'],
            initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
        (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
        start=time.monotonic();run=subprocess.run([str(executable),'input='+str((directory/'input.in').resolve())],
            env=environment,capture_output=True,text=True,timeout=max(.1,gates['maximum_total_solver_seconds']-spent))
        elapsed=time.monotonic()-start;spent+=elapsed;(directory/'run.log').write_text(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(run.stderr[-2000:])
        h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True);final=np.loadtxt(directory/'native/final.txt')
        results[label]=dict(elapsed_seconds=elapsed,initial_wall=wall_geometry(initial,512.),final_wall=wall_geometry(final,512.),
            relative_total_energy_drift=float(np.max(abs(h['total']-h['total'][0]))/h['total'][0]),
            final_axion_energy=float(h['axion_energy'][-1]))
        print(label,'completed',elapsed,'seconds',flush=True)
        if label=='empty-wall':
            results[label]['zero_axion_exact']=bool(np.all(final[:,[1,3]]==0.) and np.all(h['axion_energy']==0.) and np.all(h['window_exchange']==0.))
            continue
        t=h['t'];ea0=h['axion_energy'][0]
        integrate=lambda y:np.r_[0.,np.cumsum(.5*(y[1:]+y[:-1])*np.diff(t))]
        flux=h['window_axion_flux'];q=h['window_exchange'];out=integrate(np.maximum(flux,0.));returned=integrate(np.maximum(-flux,0.))
        residual=h['window_axion_energy']-h['window_axion_energy'][0]+integrate(flux)-integrate(q)
        indices=np.arange(0,len(t),2)
        if indices[-1]!=len(t)-1:indices=np.r_[indices,len(t)-1]
        cadence=max(abs(np.trapezoid(v,t)-np.trapezoid(v[indices],t[indices])) for v in (flux,q))/ea0
        p=Parameters(dv=0.);phi,a,pp,pa=final.T;va=terms(phi,a,p)[1];dx=.25
        diff=np.roll(a,-1)-a;density=.5*pa**2+va+.25*(diff**2+np.roll(diff,1)**2)/dx**2
        results[label].update(initial_axion_energy=float(ea0),outward_crossing_fraction=float(out[-1]/ea0),
            returned_crossing_fraction=float(returned[-1]/ea0),net_axion_exchange_fraction=float(integrate(q)[-1]/ea0),
            max_axion_continuity_relative=float(np.max(abs(residual))/ea0),output_cadence_relative=float(cadence),
            final_axion_total_ratio=float(h['axion_energy'][-1]/ea0),
            final_window_ratio=float(h['window_axion_energy'][-1]/ea0),
            final_true_vacuum_axion_fraction=float(density[phi>.9].sum()*dx/ea0),
            final_wall_region_axion_fraction=float(density[(phi>=.1)&(phi<=.9)].sum()*dx/ea0),
            final_false_vacuum_axion_fraction=float(density[phi<.1].sum()*dx/ea0))
    r=results['packet'];checks=dict(empty_exact=results['empty-wall']['zero_axion_exact'],
        axion_continuity=r['max_axion_continuity_relative']<gates['axion_continuity_relative'],
        output_cadence=r['output_cadence_relative']<gates['output_cadence_relative'],
        incident_crossing=r['outward_crossing_fraction']>gates['incident_crossing_fraction_min'],
        final_packet_coverage=abs(r['final_axion_total_ratio']-r['final_window_ratio'])<gates['final_window_missing_fraction_max'])
    report=dict(kind=__doc__,contract=contract,results=results,total_solver_seconds=spent,
        checks=checks,diagnostic_gates_passed=all(checks.values()),
        physical_scattering_validated=False,remaining=['Separate time and spatial refinement','Half-amplitude test',
            'Boundary/empty-wall radiation controls','Finite-band reference comparison','Account for residual packet energy outside detector'],
        build=build,source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('run_planar_packet_pilot.py','paper_reference.py','planar_diagnostics.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(results=results,total_solver_seconds=spent,checks=checks),indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
