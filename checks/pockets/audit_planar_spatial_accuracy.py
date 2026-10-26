"""Cheap stencil-consistency planning audit; no compiled solver evolution."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from planar_initial import PlanarSeed,prepare
from paper_reference import Parameters,head_on_reflection
from planar_stencils import laplacian,symbol


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    gates=dict(min_second_order_ratio=3.,min_fourth_order_ratio=12.,minimum_midgrid_residual_improvement=10.)
    (args.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    p=Parameters(dv=0.);rows=[]
    for speed in (0.,.3):
        for n in (2048,4096,8192):
            x,state=prepare(n,PlanarSeed(speed=speed,amplitude=0.),p);phi=state[:,0];dx=512./n
            u_phi=p.lam/2*phi*(1-phi)*(1-2*phi)
            rows.append(dict(speed=speed,n=n,dx=dx,second_residual=float(np.max(abs((1-speed**2)*laplacian(phi,dx,2)-u_phi))),
                fourth_residual=float(np.max(abs((1-speed**2)*laplacian(phi,dx,4)-u_phi)))))
    checks={}
    for speed in (0.,.3):
        r=[row for row in rows if row['speed']==speed]
        checks[str(speed)+'_second_order']=all(r[i]['second_residual']/r[i+1]['second_residual']>gates['min_second_order_ratio'] for i in (0,1))
        checks[str(speed)+'_fourth_order']=all(r[i]['fourth_residual']/r[i+1]['fourth_residual']>gates['min_fourth_order_ratio'] for i in (0,1))
        checks[str(speed)+'_improvement']=r[1]['second_residual']/r[1]['fourth_residual']>gates['minimum_midgrid_residual_improvement']
    reflected=head_on_reflection(.5,.1,.3)['reflected_lab_frequency'];dispersion=[]
    for omega in (.5,reflected):
        k=np.sqrt(omega**2-.1**2)
        for dx in (.25,.125,.0625):
            dispersion.append(dict(continuum_omega=omega,k=float(k),dx=dx,
                second_omega=float(np.sqrt(.1**2+symbol(k,dx,2))),fourth_omega=float(np.sqrt(.1**2+symbol(k,dx,4)))))
    report=dict(kind=__doc__,compiled_evolution_changed=False,gates=gates,rows=rows,dispersion=dispersion,
        checks=checks,declared_checks_passed=all(checks.values()),
        limitations='A continuum sampled-kink residual is not an evolved wall. Order-4 Hamiltonian and long-bond currents must be integrated together; negative individual bond energy terms are not a pointwise physical density.',
        next='Add an optional owned adapter hook and two ghost cells, then validate compiled order-4 evolution/energy/current against independent references. Preserve order-2 defaults and archived outputs.',
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('audit_planar_spatial_accuracy.py','planar_stencils.py','planar_initial.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(rows=rows,checks=checks),indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
