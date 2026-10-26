"""Direct finite-band continuum waveform comparison on refined planar runs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from planar_initial import PlanarSeed,modes
from paper_reference import Parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',type=Path,required=True)
    parser.add_argument('--detector-reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    gates=dict(finest_scaled_waveform_error_max=.01,minimum_improvement_factor=3.)
    (args.output/'contract.json').write_text(json.dumps(dict(gates=gates,physical_core=[32.,224.],time=220.,
        comparison='Left asymptotic region, at least 32 widths from target wall; fixed incident amplitude/momentum scales'),indent=2)+'\n')
    ref=json.loads(args.detector_reference.read_text());suite=json.loads((args.suite/'report.json').read_text())
    if not ref['passed'] or not suite['declared_checks_passed']:raise RuntimeError('Prerequisite gate failed')
    seed=PlanarSeed();p=Parameters(dv=0.);indices,k,omega,weights=modes(seed,p)
    if len(k)!=len(ref['components']):raise RuntimeError('Reference spectrum differs')
    np.testing.assert_allclose(k,[r['k'] for r in ref['components']],rtol=0.,atol=0.)
    np.testing.assert_allclose(omega,[r['omega'] for r in ref['components']],rtol=0.,atol=0.)
    rows=[]
    for name in ('time-half','space-mid','space-fine'):
        source=args.suite/(name+'-packet')/'native/final.txt';state=np.loadtxt(source);n=len(state)
        x=(np.arange(n)+.5)*512/n;mask=(x>=32)&(x<=224);z=x[mask]-256.
        a=np.zeros(mask.sum());pa=np.zeros_like(a)
        for kk,ww,weight,r in zip(k,omega,weights,ref['components']):
            amplitude=complex(r['r_real'],r['r_imag'])
            wave=np.exp(1j*kk*z)+amplitude*np.exp(-1j*kk*z)
            c=seed.amplitude*weight*np.exp(1j*kk*(256-seed.centre)-1j*ww*220.)
            a+=(c*wave).real;pa+=(-1j*ww*c*wave).real
        error_a=float(np.max(abs(a-state[mask,1]))/seed.amplitude)
        error_pa=float(np.max(abs(pa-state[mask,3]))/(seed.omega*seed.amplitude))
        rows.append(dict(case=name,n=n,dt=.0125,a_scaled_error=error_a,pi_scaled_error=error_pa,
            state_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
        np.savetxt(args.output/(name+'-comparison.csv'),np.column_stack((x[mask],state[mask,1],a,state[mask,3],pa)),
            delimiter=',',header='x,compiled_a,reference_a,compiled_pi_a,reference_pi_a',comments='')
    errors=[max(r['a_scaled_error'],r['pi_scaled_error']) for r in rows]
    factors=[errors[i]/errors[i+1] for i in (0,1)]
    checks=dict(finest_error=errors[-1]<gates['finest_scaled_waveform_error_max'],
        improves_on_refinement=all(f>gates['minimum_improvement_factor'] for f in factors))
    result=dict(kind=__doc__,scope='Stationary weak-field planar packet only; nonlinear recoil is omitted by reference',
        gates=gates,rows=rows,improvement_factors=factors,checks=checks,declared_checks_passed=all(checks.values()),
        suite_report_sha256=hashlib.sha256((args.suite/'report.json').read_bytes()).hexdigest(),
        reference_report_sha256=hashlib.sha256(args.detector_reference.read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('check_planar_reference_convergence.py','planar_initial.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
