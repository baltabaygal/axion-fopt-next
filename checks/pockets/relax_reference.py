"""Independent radial, time-averaged initial-wall preparation reference.

Minimises wall/vacuum energy plus sum_j (2l_j+1) I_j omega_j[phi]
at fixed classical mode actions I_j. This is a spherical ensemble reference,
not a lattice simulation, a thermalisation model, or a stability proof.
The potential is the unchanged phi^2 coupling used by the compiled backend.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize, root

from model import dilaton_potential, force, radial_modes


class RadialEnergy:
    def __init__(self, h, rmax, groups, lam=2., dv=.01, mass=2.):
        self.h=h; self.r=h*np.arange(1,round(rmax/h)); self.groups=groups
        self.lam=lam;self.dv=dv;self.mass=mass
        self.volume=4*np.pi*self.r**2*h
        self.faces=h*(np.arange(len(self.r)+1)+.5)
        self.off=-np.ones(len(self.r)-1)/h**2

    def evaluate(self, phi, return_modes=False):
        h,r=self.h,self.r
        padded=np.concatenate(([phi[0]],phi,[1.]))
        delta=np.diff(padded)
        wall=float(np.sum(4*np.pi*self.faces**2*delta**2/(2*h))
                   +np.sum(self.volume*dilaton_potential(phi,self.lam,self.dv)))
        flux=4*np.pi*self.faces**2*delta/h
        gradient=flux[:-1]-flux[1:]+self.volume*force(phi,0.,self.lam,self.dv,0.)[0]
        gas=0.; modes=[]
        for l, entries in self.groups.items():
            maximum=max(e['n'] for e in entries)
            diag=2/h**2+l*(l+1)/r**2+self.mass**2*phi**2
            eigenvalues,vectors=eigh_tridiagonal(diag,self.off,select='i',select_range=(0,maximum))
            for entry in entries:
                j=entry['n']; omega=np.sqrt(eigenvalues[j]); v=vectors[:,j]
                action=entry['action']; multiplicity=2*l+1
                gas+=multiplicity*action*omega
                gradient+=multiplicity*action*self.mass**2*phi*v**2/omega
                if return_modes:
                    modes.append(dict(l=l,n=j,omega=float(omega),action=action,u=v/np.sqrt(h)))
        if return_modes:
            return wall+gas,gradient,dict(wall=wall,gas=float(gas),modes=modes)
        return wall+gas,gradient


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,default=Path('outputs/pockets/balanced-pilot-01'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--h',type=float,default=.1)
    args=p.parse_args()
    if not .02<=args.h<=.2:
        p.error('h must be .02..0.2')
    args.output.mkdir(parents=True,exist_ok=False)
    baseline=json.loads((args.baseline/'contract.json').read_text())
    modes=radial_modes(20.,2.,2.,.15,.6,31.,h=.01)
    count=sum(2*m['l']+1 for m in modes)
    energy_per_mode=baseline['target_axion_energy']/count
    groups={}
    for mode in modes:
        groups.setdefault(mode['l'],[]).append(dict(n=mode['n'],action=energy_per_mode/mode['omega']))
    system=RadialEnergy(args.h,40.,groups)
    r=system.r; phi0=.5*(1+np.tanh((r-20)/2))
    (args.output/'contract.json').write_text(json.dumps(dict(
        h=args.h,rmax=40.,groups=groups,initial_radius=20.,width=2.,
        potential='unchanged lambda=2,dV=.01,mTV=2,phi^2 mass',
        baseline=str(args.baseline.resolve()),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Fixed-action spherical ensemble preparation, not evolved fields or a stability claim.',
        gates=dict(gradient_relative=2e-6,residual_density_max=1e-4,bound_frequency_margin=.1)),indent=2)+'\n')
    # Check the spectral Hellmann-Feynman force independently by directional
    # finite differences before minimising. No force fitting.
    direction=np.random.default_rng(7).normal(size=len(r));direction/=np.linalg.norm(direction)
    e0,g0=system.evaluate(phi0); eps=.01
    # Richardson cancellation avoids subtracting nearly identical spectral
    # energies at tiny steps, while removing the O(eps^2) truncation term.
    slopes=[(system.evaluate(phi0+step*direction)[0]-system.evaluate(phi0-step*direction)[0])/(2*step)
            for step in (eps,eps/2)]
    numerical=(4*slopes[1]-slopes[0])/3
    exact=float(g0@direction)
    gradient_error=abs(numerical-exact)/max(abs(exact),1.)
    if gradient_error>2e-6:
        raise RuntimeError(f'Spectral force gradient check failed: {gradient_error}')
    def objective(q):
        energy,gradient=system.evaluate(1+q/r)
        return energy,gradient/r
    result=minimize(objective,r*(phi0-1),jac=True,method='L-BFGS-B',
                    bounds=list(zip(-1.05*r,.05*r)),
                    options=dict(maxiter=4000,ftol=1e-14,gtol=1e-8,maxls=40))
    phi=1+result.x/r
    # Energy minimisation can stop with a small absolute error concentrated
    # in the tiny central cells. Polish the force-density residual itself.
    before_energy=system.evaluate(phi)[0]
    polish=root(lambda f: system.evaluate(f)[1]/system.volume,phi,
                method='hybr',options=dict(xtol=1e-9,maxfev=4000))
    polished_energy,polished_gradient=system.evaluate(polish.x)
    accepted=bool(polish.success and polished_energy<=before_energy+1e-6
                  and np.all(polish.x>=-.05) and np.all(polish.x<=1.05))
    if accepted:
        phi=polish.x
    energy,gradient,details=system.evaluate(phi,return_modes=True)
    residual=gradient/system.volume
    crossing=np.flatnonzero((phi[:-1]<.5)&(phi[1:]>=.5))
    radius=float(r[crossing[-1]]+(.5-phi[crossing[-1]])/(phi[crossing[-1]+1]-phi[crossing[-1]])*args.h) if len(crossing) else None
    max_residual=float(np.max(np.abs(residual)))
    frequencies=[m['omega'] for m in details['modes']]
    np.savez(args.output/'profile.npz',r=r,phi=phi,initial_phi=phi0,residual=residual,
             **{f"u_l{m['l']}_n{m['n']}":m['u'] for m in details['modes']})
    for m in details['modes']:
        del m['u']
    report=dict(optimizer_success=bool(result.success),message=str(result.message),iterations=result.nit,
                force_polish_success=accepted,force_polish_message=str(polish.message),
                evaluations=result.nfev,initial_energy=e0,final_energy=float(energy),
                gradient_check=gradient_error,max_force_density_residual=max_residual,
                equilibrium_radius=radius,highest_frequency=max(frequencies),**details,
                reference_pass=bool(result.success and max_residual<1e-4 and max(frequencies)<1.9),
                interpretation='Time-averaged stationary reference only. A random field realisation still fluctuates.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='modes'},indent=2))


if __name__=='__main__':
    main()
