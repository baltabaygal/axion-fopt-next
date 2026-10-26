"""Compare a saved stationary pilot to a finite-band continuum field packet.

No lattice evolution. The continuum wall is fixed/degenerate; the saved
compiled pilot has fully coupled fields and a sampled wall transient.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from planar_initial import PlanarSeed,modes
from paper_reference import Parameters
from formation_scattering_reference import scatter


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.pilot/'report.json').read_text());contract=report['contract']
    n=contract['n'];box=contract['box'];time=contract['duration'];dx=box/n
    x=(np.arange(n)+.5)*dx;z=x-box/2;seed=PlanarSeed();p=Parameters(dv=0.)
    indices,k,omega,weights=modes(seed,p);a=np.zeros(n);pa=np.zeros(n);ax=np.zeros(n);components=[]
    for j,(kk,ww,weight) in enumerate(zip(k,omega,weights)):
        row=scatter(float(ww),extent=16.,samples=513)
        r=complex(row['reflection_real'],row['reflection_imag']);kr=np.sqrt(complex(ww*ww-1.))
        wave=np.empty(n,dtype=complex);derivative=np.empty(n,dtype=complex)
        left=z< -16.;right=z>16.;middle=~(left|right)
        wave[left]=np.exp(1j*kk*z[left])+r*np.exp(-1j*kk*z[left])
        derivative[left]=1j*kk*(np.exp(1j*kk*z[left])-r*np.exp(-1j*kk*z[left]))
        wave[right]=row['wave'][-1]*np.exp(1j*kr*(z[right]-16.));derivative[right]=1j*kr*wave[right]
        wave[middle]=np.interp(z[middle],row['x'],row['wave'])
        derivative[middle]=np.interp(z[middle],row['x'],row['derivative'])
        coefficient=seed.amplitude*weight*np.exp(1j*kk*(box/2-seed.centre)-1j*ww*time)
        a+=(coefficient*wave).real;pa+=(-1j*ww*coefficient*wave).real;ax+=(coefficient*derivative).real
        components.append(dict(k=float(kk),omega=float(ww),weight=float(weight),R=row['reflected_flux'],T=row['transmitted_flux'],flux_residual=row['flux_residual']))
    phi=.5*(1+np.tanh(z));mass2=p.mass_false**2*np.exp(p.kappa*phi)
    # Evaluate the same discrete endpoint-sharing energy operator on the
    # continuum fields, to separate waveform comparison from integration conventions.
    diff=np.roll(a,-1)-a;density=.5*pa*pa+.5*mass2*a*a+.25*(diff*diff+np.roll(diff,1)**2)/dx**2
    wspec=contract['window'];d=(x-wspec['centre']+box/2)%box-box/2
    s=np.clip((abs(d)-wspec['inner'])/(wspec['outer']-wspec['inner']),0,1)
    window=1-s**3*(10+s*(-15+6*s))
    final=np.loadtxt(args.pilot/'packet/native/final.txt')
    np.savetxt(args.output/'reference-final.csv',np.column_stack((x,a,pa,ax,density)),delimiter=',',header='x,a,pi_a,gradient_a,axion_density',comments='')
    initial_energy=report['results']['packet']['initial_axion_energy']
    result=dict(kind=__doc__,scope='Single coarse pilot; diagnostic coverage failure remains unchanged',
        modes=components,parameters=dict(n=n,box=box,time=time,seed=vars(seed)),
        predicted_stationary_spectral_reflection=float(sum(ww**2*c**2*r['R'] for ww,c,r in zip(omega,weights,components))/sum(omega**2*weights**2)),
        reference_final_window_ratio=float(np.sum(window*density)*dx/initial_energy),
        reference_final_total_ratio=float(np.sum(density)*dx/initial_energy),
        max_a_error_over_incident_amplitude=float(np.max(abs(final[:,1]-a))/seed.amplitude),
        max_pi_error_over_incident_scale=float(np.max(abs(final[:,3]-pa))/(seed.omega*seed.amplitude)),
        measured_final_window_ratio=report['results']['packet']['final_window_ratio'],
        original_coverage_gate_passed=report['checks']['final_packet_coverage'],
        interpretation='Finite detector coverage is not loss. Fixed continuum wall omits coupled recoil and lattice-wall transients. No convergence gate is passed by this comparison.',
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('compare_planar_packet_reference.py','formation_scattering_reference.py','planar_initial.py','paper_reference.py')},
        pilot_report_sha256=hashlib.sha256((args.pilot/'report.json').read_bytes()).hexdigest())
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('modes','source_sha256')},indent=2))

if __name__=='__main__':main()
