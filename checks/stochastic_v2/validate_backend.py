"""Independent small-grid validation before the single-bubble pilot."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import PARAMS,read_state,rhs,run,terms
from build_backend import sha

def verify_build(exe):
    manifest=exe.parent/'build.json';meta=json.loads(manifest.read_text())
    if sha(exe)!=meta['executable_sha256']:raise RuntimeError('Executable changed')
    for name,digest in {**meta['source_sha256'],**meta['library_sha256']}.items():
        if sha(name)!=digest:raise RuntimeError('Build dependency changed: '+name)
    return meta

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();build=verify_build(args.executable);out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    gates=dict(finest_error=2e-5,min_ratio=3.5,snapshot_energy=1e-11)
    (out/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    cases=[]
    for homogeneous in (True,False):
        n=4;box=6.;x=(np.arange(n)+.5)/n*2*np.pi
        y=np.zeros((n,n,n,4));y[:]=[.4,.17,.03,-.02]
        if not homogeneous:
            y[...,0]+=.07*np.cos(x)[:,None,None];y[...,1]+=.03*np.sin(x)[None,:,None]
            y[...,2]+=.02*np.cos(x)[None,None,:]
        reference=solve_ivp(lambda t,v:rhs(t,v,n,box),(0,.4),y.ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
        if not reference.success:raise RuntimeError(reference.message)
        np.save(out/('reference-'+str(homogeneous)+'.npy'),reference.y[:,-1].reshape(y.shape))
        errors=[];runs=[]
        for dt in (.02,.01,.005):
            directory=out/f'{homogeneous}-{dt}';seconds=run(args.executable,directory,y,box,dt,.4,cadence=.02,frames=.4)
            _,initial=read_state(directory/'native/initial.bin');_,final=read_state(directory/'native/final.bin')
            h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True)
            error=float(np.max(np.abs(final-reference.y[:,-1].reshape(y.shape))));errors.append(error)
            energy=sum(float(v.sum()) for v in terms(final,box/n).values())*(box/n)**3
            runs.append(dict(dt=dt,error=error,initial_exact=bool(np.array_equal(initial,y)),
                energy_error=float(abs(energy-h['total'][-1])/abs(h['total'][-1])),seconds=seconds))
        ratios=[errors[i]/errors[i+1] for i in (0,1)]
        cases.append(dict(homogeneous=homogeneous,runs=runs,ratios=ratios,
            passed=bool(errors[-1]<gates['finest_error'] and min(ratios)>gates['min_ratio'] and
            all(r['initial_exact'] and r['energy_error']<gates['snapshot_energy'] for r in runs))))
    # Global canonical Hamiltonian derivative equals minus the compiled-model
    # force reference, including periodic bond energy and full reciprocal force.
    rng=np.random.default_rng(42);y=rng.normal(size=(4,4,4,4))*.15;y[...,0]+=.4
    direction=rng.normal(size=y.shape);direction[...,2:]=0;eps=1e-6;dx=1.5
    energy=lambda state:sum(float(v.sum()) for v in terms(state,dx).values())*dx**3
    derivative=(energy(y+eps*direction)-energy(y-eps*direction))/(2*eps)
    acceleration=rhs(0,y.ravel(),4,6.).reshape(y.shape)[...,2:]
    expected=-float(np.sum(acceleration*direction[...,:2]))*dx**3
    gradient_error=abs(derivative-expected)/max(1,abs(expected))
    report=dict(cases=cases,gradient_error=gradient_error,build=build,
        passed=bool(all(c['passed'] for c in cases) and gradient_error<1e-7),
        source_sha256={str(Path(__file__).resolve()):sha(__file__),str(Path(__file__).with_name('model.py').resolve()):sha(Path(__file__).with_name('model.py'))})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(passed=report['passed'],cases=cases,gradient_error=gradient_error),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
