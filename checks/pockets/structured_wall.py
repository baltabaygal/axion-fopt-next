"""Independent preparation for the structured toy dilaton wall.

This is a radial/kink reference, not a compiled lattice evolution.
Potential: STRUCTURED_WALL_EXPERIMENT.md. The transformed force and
curvature polynomials provide checks independent of a factored C++ form.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq
from scipy.special import expit


def check_shape(b):
    if not np.isfinite(b) or not 0<=b<16:
        raise ValueError('Controlled structured wall requires 0<=b<16')


def potential(phi,lam,dv,b):
    check_shape(b)
    return (1-phi)**2*(lam/4*phi**2*(1+b*(2*phi-1)**8)+dv*(1+2*phi+3*phi**2))


def force(phi,lam,dv,b):
    check_shape(b);y=2*phi-1
    return lam*phi*(1-phi)*y*(-.5+b*y**6*(1-1.5*y*y))+12*dv*phi*phi*(phi-1)


def curvature(phi,lam,dv,b):
    check_shape(b);y=2*phi-1
    return lam/4*(-1+3*y*y+14*b*y**6-45*b*y**8+33*b*y**10)+dv*(36*phi*phi-24*phi)


def kink(lam,b,extent=30.):
    """Exact degenerate first-integral ODE in logit coordinates; stable tails."""
    check_shape(b)
    def rhs(_,s):
        return np.sqrt(lam/2)*np.sqrt(1+b*(2*expit(s)-1)**8)
    positive=solve_ivp(rhs,(0,extent),[0.],method='DOP853',dense_output=True,rtol=2e-13,atol=2e-14,max_step=.1)
    negative=solve_ivp(rhs,(0,-extent),[0.],method='DOP853',dense_output=True,rtol=2e-13,atol=2e-14,max_step=.1)
    if not positive.success or not negative.success:
        raise RuntimeError('Kink reference integration failed')
    def profile(z):
        original=np.asarray(z,dtype=float);flat=original.ravel();values=np.empty_like(flat)
        pos=(flat>=0)&(flat<=extent);neg=(flat<0)&(flat>=-extent)
        values[pos]=expit(positive.sol(flat[pos])[0]);values[neg]=expit(negative.sol(flat[neg])[0])
        values[flat>extent]=1.;values[flat<-extent]=0.
        return values.reshape(original.shape)
    return profile


def tension(lam,b):
    check_shape(b)
    return quad(lambda p:np.sqrt(lam/2)*p*(1-p)*np.sqrt(1+b*(2*p-1)**8),0,1,epsabs=1e-13,epsrel=1e-13)[0]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    lam,dv,b=2.,.01,12.
    declaration=Path(__file__).with_name('STRUCTURED_WALL_EXPERIMENT.md')
    (args.output/'contract.json').write_text(json.dumps(dict(lambda_=lam,delta_v=dv,b=b,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        declaration_sha256=hashlib.sha256(declaration.read_bytes()).hexdigest(),
        gates=dict(derivative_relative=2e-7,curvature_relative=2e-7,kink_energy_relative=1e-7,
                   kink_residual=.005,kink_refinement=3.5,interior_barriers=1)),indent=2)+'\n')
    phi=np.random.default_rng(31).uniform(-.05,1.05,100);eps=1e-6
    fd=(potential(phi+eps,lam,dv,b)-potential(phi-eps,lam,dv,b))/(2*eps)
    fh=(force(phi+eps,lam,dv,b)-force(phi-eps,lam,dv,b))/(2*eps)
    gradient_error=float(np.max(np.abs(fd-force(phi,lam,dv,b))/np.maximum(1.,np.abs(fd))))
    hessian_error=float(np.max(np.abs(fh-curvature(phi,lam,dv,b))/np.maximum(1.,np.abs(fh))))
    x=np.linspace(1e-9,1-1e-9,20001)
    barriers=[brentq(lambda q:force(q,lam,dv,b),a,c) for a,c in zip(x[:-1],x[1:])
              if force(a,lam,dv,b)*force(c,lam,dv,b)<0]
    profile=kink(lam,b);sigma=tension(lam,b);z=np.linspace(-15,15,30001);f=profile(z)
    derivative=np.sqrt(lam/2)*f*(1-f)*np.sqrt(1+b*(2*f-1)**8)
    energy=float(np.trapezoid(.5*derivative**2+potential(f,lam,0.,b),z))
    residuals=[]
    for dz in (.02,.01):
        positions=np.arange(-6,6+dz/2,dz);values=profile(positions)
        residual=(values[2:]-2*values[1:-1]+values[:-2])/dz**2-force(values[1:-1],lam,0.,b)
        residuals.append(float(np.max(np.abs(residual))))
    width=quad(lambda q:1/(np.sqrt(lam/2)*q*(1-q)*np.sqrt(1+b*(2*q-1)**8)),.1,.9,epsabs=1e-12)[0]
    report=dict(gradient_error=gradient_error,hessian_error=hessian_error,interior_barriers=barriers,
        sigma=sigma,kink_energy=energy,kink_energy_relative_error=abs(energy/sigma-1),
        kink_residuals=residuals,kink_refinement=residuals[0]/residuals[1],width_10_90=width,
        false_mass=float(np.sqrt(curvature(0.,lam,dv,b))),true_mass=float(np.sqrt(curvature(1.,lam,dv,b))),
        passed=bool(gradient_error<2e-7 and hessian_error<2e-7 and len(barriers)==1
            and abs(energy/sigma-1)<1e-7 and max(residuals)<.005 and residuals[0]/residuals[1]>3.5),
        interpretation='Independent potential/kink reference only; compiled implementation and field survival not established.')
    np.savez(args.output/'kink.npz',z=z,phi=f,derivative=derivative)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
