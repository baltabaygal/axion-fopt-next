"""Tiny standalone C++ RK4 vs independent SciPy DOP853; no lattice claims."""
import argparse
from dataclasses import asdict
import hashlib
import io
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import numpy as np
from scipy.integrate import solve_ivp
from paper_reference import Parameters,force,terms

HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    compiler=shutil.which('c++')
    if not compiler:raise RuntimeError('C++ compiler required')
    binary=args.output/'probe'
    command=[compiler,'-std=c++17','-O0',str(HERE/'formation_homogeneous_probe.cpp'),'-o',str(binary)]
    build=subprocess.run(command,capture_output=True,text=True,timeout=30)
    (args.output/'build.log').write_text(build.stdout+build.stderr)
    build.check_returncode()
    p=Parameters();times=np.linspace(0.,4.,21);y0=[0.,2.,0.,.15,0.]
    def rhs(t,y):
        f,g=force(y[0],y[1],p);_,va,_,_=terms(y[0],y[1],p)
        return [y[2],y[3],-float(f),-float(g),float(p.kappa*va*y[2])]
    ref=solve_ivp(rhs,(0.,4.),y0,method='DOP853',rtol=2e-13,atol=2e-15,t_eval=times)
    if not ref.success:raise RuntimeError(ref.message)
    np.savetxt(args.output/'reference.csv',np.column_stack((times,ref.y.T)),delimiter=',',
               header='t,phi,a,phi_dot,a_dot,integrated_Q',comments='')
    # Gates declared before comparison, using fixed initial amplitude scales.
    gates=dict(fine_scaled_state_error=2e-7,fine_relative_energy_drift=2e-7,
               fine_relative_exchange_residual=2e-7,min_error_reduction=10.)
    results=[]
    scales=np.array([1.,2.,1.,.2,1.])
    for steps in (80,160,320):
        run=subprocess.run([str(binary.resolve()),str(steps)],capture_output=True,text=True,check=True,timeout=10)
        (args.output/f'rk4-{steps}.csv').write_text(run.stdout)
        v=np.loadtxt(io.StringIO(run.stdout),delimiter=',',skiprows=1)
        np.testing.assert_allclose(v[:,0],times,rtol=0.,atol=1e-14)
        results.append(dict(steps=steps,dt=4./steps,
            scaled_state_error=float(np.max(np.abs(v[:,1:6]-ref.y.T)/scales)),
            relative_energy_drift=float(np.max(np.abs(v[:,6]-v[0,6]))/v[0,6]),
            relative_exchange_residual=float(np.max(np.abs(v[:,7]-v[0,7]-v[:,5]))/v[0,7])))
    fine=results[-1]
    checks={name:fine[name.removeprefix('fine_')]<limit for name,limit in gates.items() if name.startswith('fine_')}
    checks['fourth_order_time_refinement']=all(results[i]['scaled_state_error']/results[i+1]['scaled_state_error']>gates['min_error_reduction'] for i in (0,1))
    checks['false_vacuum_backreaction']=bool(np.min(ref.y[0])<-.001)
    sources=['formation_potential.h','formation_homogeneous_probe.cpp','paper_reference.py','validate_formation_homogeneous.py']
    report=dict(kind='standalone homogeneous kernel validation; not lattice evolution',
        parameters=asdict(p),initial_state=y0,time='Minkowski cosmic time; t=0..4',
        reference=dict(method='SciPy DOP853',rtol=2e-13,atol=2e-15,nfev=ref.nfev),
        normalization=dict(state_scales=scales.tolist(),energy='initial total or axion energy as appropriate'),
        gates=gates,results=results,checks=checks,passed=all(checks.values()),
        provenance=dict(python=sys.version,platform=platform.platform(),build_command=command,
            source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in sources},
            executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()))
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],results=results,checks=checks),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
