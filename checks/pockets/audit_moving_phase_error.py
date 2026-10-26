"""Read-only phase diagnosis using independently evolved empty-wall positions.

The effective-speed reference is a diagnostic counterfactual, not a
replacement for the declared fixed-speed continuum convergence gate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from planar_diagnostics import wall_geometry
from moving_packet_reference import left_asymptotic_fields
from planar_initial import PlanarSeed


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    suite=json.loads((args.suite/'report.json').read_text());ref=json.loads(args.reference.read_text())
    predictions=ref['mode_predictions'];rest=np.array([r['omega'] for r in predictions])
    real=np.array([r['reflection_real'] for r in predictions]);imag=np.array([r['reflection_imag'] for r in predictions])
    rows=[]
    for case,root in (('space-coarse',args.suite/'space-coarse'),('base',args.baseline),('space-fine',args.suite/'space-fine')):
        empty_path=root/'empty/native/final.txt';state_path=root/'packet/native/final.txt'
        empty=np.loadtxt(empty_path);state=np.loadtxt(state_path);n=len(state)
        geometry=wall_geometry(empty,512.);u=(256.-geometry['position'])/180.
        drift=geometry['position']-(256.-.3*180.)
        x=(np.arange(n)+.5)*512./n;core=(x>=32.)&(x<=184.)
        gamma=1/np.sqrt(1-u*u);diagnostic=[]
        for r in predictions:
            k=r['k'];om=r['lab_omega'];om_rest=gamma*(om+u*k)
            changed=dict(r,reflection_real=float(np.interp(om_rest,rest,real)),
                reflection_imag=float(np.interp(om_rest,rest,imag)),
                reflected_lab_omega=float(gamma**2*(om*(1+u*u)+2*u*k)))
            diagnostic.append(changed)
        seed=PlanarSeed(speed=u)
        a,pa=left_asymptotic_fields(x[core],180.,diagnostic,seed)
        counterfactual=max(float(np.max(abs(state[core,1]-a))/.001),float(np.max(abs(state[core,3]-pa))/.0005))
        original=suite['cases'][case]['scaled_reference_field_error']
        rows.append(dict(case=case,n=n,dx=512./n,empty_wall_position=geometry['position'],
            empty_wall_drift_from_continuum=drift,empty_average_approach_speed=u,
            original_fixed_speed_scaled_error=original,effective_speed_diagnostic_scaled_error=counterfactual,
            error_reduction_factor=original/counterfactual,
            state_sha256=hashlib.sha256(state_path.read_bytes()).hexdigest(),empty_sha256=hashlib.sha256(empty_path.read_bytes()).hexdigest()))
    result=dict(kind=__doc__,rows=rows,solver_launched=False,
        original_direct_waveform_gate_passed=suite['checks']['direct_waveform'],
        interpretation='Empty-wall drift decreases with dx. Counterfactual uses average endpoint speed, not a measured full trajectory. Wave propagation dispersion and transients remain. Original fixed-u gate is unchanged.',
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.suite/'report.json',args.reference)},
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('audit_moving_phase_error.py','moving_packet_reference.py','planar_diagnostics.py','planar_initial.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
