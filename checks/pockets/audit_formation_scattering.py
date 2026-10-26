"""Light continuum ODE reference for initial and repeated wall encounters."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from paper_reference import head_on_reflection
from formation_scattering_reference import scatter


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    gates=dict(flux_balance=2e-8,domain_reflection_amplitude=2e-7,
               tolerance_reflection_amplitude=2e-8,subthreshold_transmitted_flux=1e-14)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    first=[head_on_reflection(omega,.1,.3) for omega in (.3,.5)]
    second=head_on_reflection(first[1]['reflected_lab_frequency'],.1,.3)
    frequencies=[first[0]['wall_frame_frequency'],first[1]['wall_frame_frequency'],
                 1.01,second['wall_frame_frequency'],1.5]
    results=[]
    for j,omega in enumerate(frequencies):
        rows=[]
        for label,extent,rtol,atol in (('base',12.,2e-11,2e-13),('domain',16.,2e-11,2e-13),
                                      ('tight',16.,2e-13,2e-15)):
            row=scatter(omega,extent=extent,rtol=rtol,atol=atol)
            curve=np.column_stack((row.pop('x'),row['wave'].real,row['wave'].imag,
                                   row['derivative'].real,row['derivative'].imag))
            row.pop('wave');row.pop('derivative')
            np.savetxt(args.output/f'omega-{j}-{label}.csv',curve,delimiter=',',
                       header='x,psi_real,psi_imag,derivative_real,derivative_imag',comments='')
            rows.append(row)
        def amplitude(r):return complex(r['reflection_real'],r['reflection_imag'])
        domain_error=abs(amplitude(rows[0])-amplitude(rows[1]))
        tolerance_error=abs(amplitude(rows[1])-amplitude(rows[2]))
        results.append(dict(omega=omega,controls=rows,domain_amplitude_error=domain_error,
            tolerance_amplitude_error=tolerance_error,
            passed=all(r['flux_residual']<gates['flux_balance'] for r in rows) and
                domain_error<gates['domain_reflection_amplitude'] and tolerance_error<gates['tolerance_reflection_amplitude'] and
                (omega>=1. or rows[-1]['transmitted_flux']<gates['subthreshold_transmitted_flux'])))
    report=dict(kind=__doc__,scope='Continuum linear weak-axion ODE only; no coupled lattice evolution or survival evidence',
        background=dict(phi='(1+tanh(x/w))/2',width=1.,lambda_=8.,delta_v=0.,
            mass_false=.1,mass_true=1.,kappa=float(np.log(100.)),
            kinetic='canonical',time='Minkowski cosmic time; psi(x) exp(-i omega t)',
            interpretation='Exact degenerate dilaton kink at zero axion; weak axion backreaction is higher order'),
        normalization='Unit incident field amplitude; flux fractions are in the stationary wall frame',
        initial_head_on_kinematics=first,second_upper_encounter=second,
        gates=gates,results=results,passed=all(r['passed'] for r in results),
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('audit_formation_scattering.py','formation_scattering_reference.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],fluxes=[dict(omega=r['omega'],R=r['controls'][-1]['reflected_flux'],
        T=r['controls'][-1]['transmitted_flux'],domain_error=r['domain_amplitude_error'],
        tolerance_error=r['tolerance_amplitude_error']) for r in results]),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
