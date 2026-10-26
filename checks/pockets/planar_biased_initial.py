"""Explicit approximate biased-wall seed and independent initial force audit.

The degenerate kink is only an initial profile when deltaV>0. No claim of
an exact traveling solution is made; planar_initial's exact-kink guard stays.
"""
from dataclasses import replace,asdict
import argparse,json
from pathlib import Path
import numpy as np
from planar_initial import prepare,PlanarSeed
from paper_reference import Parameters,force,terms
from planar_stencils import laplacian,gradient_density


def prepare_biased(n,seed,p):
    if p.dv<=0:raise ValueError('Require strictly positive vacuum bias')
    return prepare(n,seed,replace(p,dv=0.))


def audit(n=4096,seed=PlanarSeed(speed=.3,amplitude=.16),p=Parameters(dv=.0002)):
    x,s=prepare_biased(n,seed,p);dx=seed.box/n
    phi,a,pp,pa=s.T;gamma=1/np.sqrt(1-seed.speed**2)
    # Analytic derivative of the rising target profile; exponentially remote
    # partner contributes below roundoff on the integration interval.
    z=gamma*(x-seed.box/2)/seed.width
    tangent=.5*gamma/seed.width*(1-np.tanh(z)**2)
    core=abs(x-seed.box/2)<16*seed.width
    f,fa=force(phi,a,p);f0,_=force(phi,a,replace(p,dv=0.))
    bias_pressure=float(np.sum((f-f0)[core]*tangent[core])*dx)
    # P_wall=-integral(pi_phi phi_x); dP/dt=+integral(V_phi phi_x).
    # Thus target force=-deltaV for the rising FV->TV interface.
    target_force=bias_pressure
    sigma=np.sqrt(p.lam/2)/6
    expected=-p.dv
    u,v,_,_=terms(phi,a,p)
    return dict(parameters=asdict(p),seed=asdict(seed),n=n,dx=dx,
        time_coordinate='Minkowski cosmic time; scale factor 1',
        initial_profile='Boosted degenerate kink sampled as approximate biased initial data',
        target_bias_force=target_force,expected_target_bias_force=expected,
        relative_force_error=abs(target_force-expected)/p.dv,
        degenerate_rest_tension=sigma,
        thin_wall_initial_acceleration=-p.dv/(sigma*gamma**3),
        initial_acceleration_field_max=float(np.max(abs(laplacian(phi,dx,4)-f))),
        initial_axion_energy_per_area=float(np.sum(.5*pa**2+v+gradient_density(a,dx,4))*dx),
        initial_false_volume_bias_per_area=p.dv*seed.box/2,
        limitations='Force projection checks sign/normalization only; no evolution, equilibrium or lifetime evidence')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rows=[audit(n) for n in (2048,4096,8192)]
    report=dict(rows=rows,force_gate_relative=1e-8,
        checks_passed=all(r['relative_force_error']<1e-8 for r in rows),solver_runs=0)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not report['checks_passed']:raise SystemExit(1)

if __name__=='__main__':main()
