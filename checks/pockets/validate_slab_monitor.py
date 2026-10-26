"""Validate native two-interface monitor: independent endpoints, continuity, passivity."""
import argparse, json, os, subprocess, hashlib
from pathlib import Path
import numpy as np
from planar_distributed_initial import SlabSeed, prepare
from slab_retention import weights
from planar_diagnostics import wall_geometry
from planar_stencils import window_budget
from paper_reference import Parameters

HERE=Path(__file__).resolve().parent

def integral(y,t):return np.r_[0.,np.cumsum(.5*(y[1:]+y[:-1])*np.diff(t))]

def execute(exe, directory, state, box, dt, duration, padding, order=4, cadence=None):
    directory.mkdir(parents=True,exist_ok=False)
    np.savetxt(directory/'input.txt',state,fmt='%.17g')
    params=dict(N=len(state),lSide=box,dt=dt,tMax=duration,expansion=False,evolver='PV2',PS_type=2,
        line_spatial_order=order,line_slab_padding=padding,**{'lambda':8.},delta_v=.0002,
        mass_fv=.1,mass_kappa=float(np.log(100)),axion_fa=10.,output_every=round((cadence or dt)/dt),
        initial_state=str((directory/'input.txt').resolve()),result_dir=str((directory/'native').resolve()))
    (directory/'input.in').write_text(''.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}\n' for k,v in params.items()))
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    run=subprocess.run([str(exe),'input='+str((directory/'input.in').resolve())],env=env,capture_output=True,text=True,timeout=30)
    (directory/'run.log').write_text(run.stdout+run.stderr)
    if run.returncode:raise RuntimeError(run.stderr[-2000:])
    return np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True),np.loadtxt(directory/'native/final.txt')

def reference(s,box,p,padding,order):
    right=wall_geometry(s,box);mirror=wall_geometry(s[::-1],box)
    left=box-mirror['position'];lv=-mirror['phase_velocity'];rv=right['phase_velocity']
    x=(np.arange(len(s))+.5)*box/len(s)
    w,wd=weights(x,box,(left+right['position'])/2,(right['position']-left)/2,(lv+rv)/2,(rv-lv)/2,padding)
    return window_budget(s,box,w,wd,p,order)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--executable',type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False);exe=a.executable.resolve()
    gates=dict(reference_absolute=1e-11,finest_balance=1e-5,minimum_ratio=3.,passive_exact=True)
    (a.output/'contract.json').write_text(json.dumps(gates,indent=2)+'\n')
    p=Parameters(dv=.0002);n=512;box=128.
    x,s,_=prepare(n,SlabSeed(box=box,length=64.,taper=8.))
    # Nonzero field at both interfaces exercises flux and exchange immediately.
    s[:,2]-=.02*(np.roll(s[:,0],-1)-np.roll(s[:,0],1))/(2*box/n)
    s[:,1]=.02*np.sin(2*np.pi*x/box);s[:,3]=.003*np.cos(4*np.pi*x/box)
    rows=[]
    for order in (2,4):
        for dt in (.01,.005,.0025):
            d=a.output/f'order-{order}-dt-{dt}';h,f=execute(exe,d,s,box,dt,.2,4.,order)
            agreement=max(abs(v-h['slab_'+k][i]) for i,state in ((0,s),(-1,f)) for k,v in reference(state,box,p,4.,order).items())
            t=h['t'];q=integral(h['slab_exchange'],t)
            residual=max(float(np.max(abs(h['slab_'+sector+'_energy']-h['slab_'+sector+'_energy'][0]+integral(h['slab_'+sector+'_flux'],t)+sign*q))) for sector,sign in (('phi',1),('axion',-1)))/h['slab_energy'][0]
            off,offstate=execute(exe,a.output/f'off-{order}-{dt}',s,box,dt,.2,0.,order)
            passive=bool(np.array_equal(f,offstate) and all(np.array_equal(h[k],off[k]) for k in off.dtype.names))
            rows.append(dict(order=order,dt=dt,reference_absolute=float(agreement),balance=residual,passive_exact=passive))
    checks=dict(reference=all(r['reference_absolute']<gates['reference_absolute'] for r in rows),passive=all(r['passive_exact'] for r in rows),
        continuity=all(rows[i]['balance']<gates['finest_balance'] for i in (2,5)),convergence=all(rows[i]['balance']/rows[i+1]['balance']>3 for i in (0,1,3,4)))
    result=dict(gates=gates,rows=rows,checks=checks,passed=all(checks.values()),build=json.loads((exe.parent/'build.json').read_text()),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (a.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(rows=rows,checks=checks),indent=2))
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
