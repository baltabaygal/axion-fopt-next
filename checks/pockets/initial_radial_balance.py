"""Time-averaged radial force for the actual seeded classical mode actions.

A spherical spectral diagnostic, not a field evolution or local isotropy
assumption. The tanh profile is varied in radius at fixed width/actions.
Only the unmodified vacuum potential is supported here.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize_scalar
from model import radial_modes,dilaton_potential,force,lambda_for_width


def analyse(state_contract):
    spec=state_contract['spec']
    if spec.get('vacuum_stiffening',0.):raise ValueError('Tanh reference requires b=0')
    R,w,h=spec['r0'],spec['width'],spec['radial_h']
    mass,power=spec['mass_tv'],spec['mass_power'];lam=lambda_for_width(w);dv=spec['delta_v']
    modes=radial_modes(R,w,mass,spec['omega_min'],spec['omega_max'],spec['radial_rmax'],h=h,mass_power=power)
    if not modes:raise ValueError('No trapped modes')
    recorded=state_contract['modes']
    if len(recorded)!=len(modes) or any((a['l'],a['n'])!=(b['l'],b['n']) or abs(a['omega']-b['omega'])>1e-10 for a,b in zip(recorded,modes)):
        raise ValueError('Reference modes do not match archived initial state')
    rng=np.random.default_rng(spec['seed']);amplitude=state_contract['amplitude']
    actions=[]
    # Match eigen_gas's group, radial-mode, real-harmonic random draw order.
    for l in sorted({m['l'] for m in modes}):
        for m in (m for m in modes if m['l']==l):
            energy=0.
            for _ in range(2*l+1):
                q,p=rng.normal(),rng.normal();energy+=.5*amplitude**2*(q*q+p*p)
            actions.append((m,energy/m['omega']))
    r=modes[0]['r'];phi=.5*(1+np.tanh((r-R)/w));drphi=2/w*phi*(1-phi)
    drrphi=4/w**2*phi*(1-phi)*(1-2*phi)
    wall_derivative=float(np.trapezoid(4*np.pi*r*r*(-drphi*drrphi-force(phi,0.,lam,dv,0.)[0]*drphi),r))
    gas_derivative=0.;gas_energy=0.
    for m,action in actions:
        domega=float(np.sum(mass**2*power*phi**(power-1)*(-drphi)*m['u']**2)*h/(2*m['omega']))
        gas_derivative+=action*domega;gas_energy+=action*m['omega']
    def energy(radius,check_labels=False):
        f=.5*(1+np.tanh((r-radius)/w));df=2/w*f*(1-f)
        wall=float(np.trapezoid(4*np.pi*r*r*(.5*df*df+dilaton_potential(f,lam,dv)),r))
        gas=0.;off=-np.ones(len(r)-1)/h**2
        for l in sorted({m['l'] for m,_ in actions}):
            entries=[(m,act) for m,act in actions if m['l']==l]
            values=eigh_tridiagonal(2/h**2+l*(l+1)/r**2+mass**2*f**power,off,eigvals_only=True,select='i',select_range=(0,max(m['n'] for m,_ in entries)))
            if check_labels and any(abs(np.sqrt(values[m['n']])-m['omega'])>1e-8 for m,_ in entries):
                raise ValueError('Band excludes lower radial modes; selected-index labels are unsupported')
            gas+=sum(act*np.sqrt(values[m['n']]) for m,act in entries)
        return np.array([wall,gas])
    energy(R,check_labels=True)
    slopes=[(energy(R+eps)-energy(R-eps))/(2*eps) for eps in (.02,.01)]
    finite=(4*slopes[1]-slopes[0])/3;analytic=np.array([wall_derivative,gas_derivative])
    errors=np.abs(finite-analytic)/np.maximum(1.,np.abs(analytic))
    equilibrium=minimize_scalar(lambda radius:energy(radius).sum(),bounds=(.8*R,1.2*R),
                                method='bounded',options=dict(xatol=1e-5))
    eq_radius=float(equilibrium.x)
    eq_derivative=float((energy(eq_radius+.01).sum()-energy(eq_radius-.01).sum())/.02)
    return dict(wall_energy_derivative=wall_derivative,gas_energy_derivative=float(gas_derivative),
                net_outward_radial_force=float(-analytic.sum()),
                force_imbalance_fraction=float(-analytic.sum()/wall_derivative),
                continuum_gas_energy=float(gas_energy),lattice_gas_energy=spec['target_axion_energy'],
                dirichlet_gas_force_estimate=float(gas_energy/R),
                derivative_check_relative=errors.tolist(),derivative_gate=1e-5,
                fixed_profile_equilibrium_radius=eq_radius,
                fixed_profile_equilibrium_energy=float(equilibrium.fun),
                equilibrium_force_residual_relative=abs(eq_derivative)/wall_derivative,
                equilibrium_in_search_interior=bool(.801*R<eq_radius<1.199*R),
                passed=bool(np.max(errors)<1e-5 and equilibrium.success and abs(eq_derivative)/wall_derivative<1e-4 and .801*R<eq_radius<1.199*R),
                interpretation='Seeded fixed-action spherical mean radial force, not instantaneous local stress, stability, thermalisation or evolved equilibrium.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--state-contract',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'contract.json').write_text(json.dumps(dict(source=str(args.state_contract.resolve()),source_sha256=hashlib.sha256(args.state_contract.read_bytes()).hexdigest(),analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),derivative_gate=1e-5),indent=2)+'\n')
    report=analyse(json.loads(args.state_contract.read_text()))
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if not report['passed']:raise RuntimeError('Independent derivative check failed')


if __name__=='__main__':main()
