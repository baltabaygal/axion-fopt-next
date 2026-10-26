"""Prepare finite-band moving-wall predictions and detector checks, no solver."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from planar_initial import PlanarSeed,prepare
from paper_reference import Parameters,head_on_reflection
from moving_packet_reference import reference_modes,packet_fields


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    contract=dict(n=4096,box=512.,dt=.0125,time_end=180.,speed=.3,
        detector=dict(centre=88.,inner=72.,outer=104.),
        gates=dict(initial_coverage_min=.999,final_reflected_coverage_min=.98,minimum_window_at_encounter_max=.01,initial_energy_error_max=1e-6,initial_core_field_error_max=1e-6))
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    seed=PlanarSeed(speed=.3);p=Parameters(dv=0.);rows=reference_modes(seed,p,quadrature_factor=8)
    n=contract['n'];x=(np.arange(n)+.5)*512/n;dx=512/n
    incident=rows[0]['quadrature_period']/2*seed.amplitude**2*sum(r['lab_omega']**2*r['weight']**2 for r in rows)
    reflected=rows[0]['quadrature_period']/2*seed.amplitude**2*sum(r['lab_omega']*r['reflected_lab_omega']*r['weight']**2*r['reflected_flux'] for r in rows)
    energy_weights=np.array([r['lab_omega']*r['reflected_lab_omega']*r['weight']**2*r['reflected_flux'] for r in rows])
    frequency_energy_mean=float(np.sum(energy_weights*np.array([r['reflected_lab_omega'] for r in rows]))/energy_weights.sum())
    spec=contract['detector'];d=(x-spec['centre']+256)%512-256;z=np.clip((abs(d)-spec['inner'])/(spec['outer']-spec['inner']),0,1)
    window=1-z**3*(10+z*(-15+6*z));samples=[]
    for t in (0.,90.,95.,100.,105.,110.,180.):
        a,pa,gradient,phi=packet_fields(x,t,rows,seed)
        mass2=p.mass_false**2*np.exp(p.kappa*phi)
        density=.5*(pa*pa+gradient*gradient+mass2*a*a)
        samples.append(dict(time=t,window_over_incident=float(np.sum(window*density)*dx/incident),
            total_over_incident=float(density.sum()*dx/incident)))
        if t==0.:
            _,original=prepare(n,seed,p);core=(x>=32)&(x<=224)
            initial_core_error=float(max(np.max(abs(a[core]-original[core,1]))/seed.amplitude,np.max(abs(pa[core]-original[core,3]))/(seed.omega*seed.amplitude)))
        if t in (0.,180.):np.savetxt(args.output/f'fields-{t}.csv',np.column_stack((x,a,pa,gradient,phi)),delimiter=',',header='x,a,pi_a,gradient_a,phi',comments='')
    initial=samples[0]['window_over_incident'];final=samples[-1]['window_over_incident']/(reflected/incident)
    minimum=min(s['window_over_incident'] for s in samples[1:-1]);g=contract['gates']
    checks=dict(initial=initial>g['initial_coverage_min'],final=final>g['final_reflected_coverage_min'],
        incident_exit=minimum<g['minimum_window_at_encounter_max'],
        initial_energy=abs(samples[0]['total_over_incident']-1)<g['initial_energy_error_max'],
        initial_core=initial_core_error<g['initial_core_field_error_max'])
    result=dict(kind=__doc__,contract=contract,initial_reference_energy=incident,
        reflected_lab_energy_ratio=reflected/incident,reflected_energy_weighted_frequency=frequency_energy_mean,
        central_mode_kinematics=head_on_reflection(.5,.1,.3),samples=samples,
        quadrature_factor=8,initial_core_scaled_error=initial_core_error,final_reflected_coverage=final,checks=checks,reference_design_passed=all(checks.values()),
        normalization='Reflected lab energy uses omega_out/omega_in times wall-frame R; not unboosted R alone',
        scope='Fixed-speed degenerate weak-field reference, no coupled simulation or survival conclusion',
        mode_predictions=[{k:r[k] for k in ('k','lab_omega','omega','rest_k','weight','reflected_flux','transmitted_flux','reflected_lab_omega','reflection_real','reflection_imag')} for r in rows],
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('audit_moving_packet_reference.py','moving_packet_reference.py','formation_scattering_reference.py','planar_initial.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('mode_predictions','source_sha256')},indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
