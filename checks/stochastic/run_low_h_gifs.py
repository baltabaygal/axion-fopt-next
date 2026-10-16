"""Representative H*=0.05 animation run with mt=1 and fixed Hubble geometry."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from scipy.integrate import solve_ivp
from run_field import ROOT,read_state
from analytic import Case,survival
from axion_fopt_next import Catalogue
from axion_fopt_next.interchange import write_arrival_map


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/stochastic/representative-H0.05')
    args=parser.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    h=.05;box=60.;width=3.;origin=20.;end=60.
    cfg=dict(h=h,beta=h,box=box,width=width,origin=origin,n=64,birth_end=50.,field_end=end,dt=.01,mass_true=1.,mass_false=0.,speed=.5)
    (root/'animation.json').write_text(json.dumps(cfg,indent=2)+'\n')
    z=np.load(ROOT/'outputs/stochastic/field-reference/catalogue.npz')
    catalogue=Catalogue(box,.5,z['times']/h,z['positions']/h)
    np.savez_compressed(root/'catalogue.npz',times=catalogue.times,positions=catalogue.positions)
    old=Case('old',3.,.5,2.5);case=Case('lowH',box,.5,50.,h=h,beta=h)
    np.testing.assert_allclose(survival(case,np.linspace(0,50,26)),survival(old,np.linspace(0,2.5,26)),atol=1e-12,rtol=0)
    (root/'EXPERIMENT.md').write_text('''# H*=0.05 representative animation contract

Keep mt=omega_star=f=1,mf=0 and beta/H=1,H*L=3,H*width=0.15.
Thus beta=0.05,L=60,width=3. Rescale the validated seed-0 614-bubble
catalogue by 1/H (positions and conformal births), preserving its normalized
geometry. This does not reuse/rescale the old axion field: solve it anew.
Birth window tau=0..50, theta=1,theta'=0 imposed at tau=20,t=10,a=1;
evolve to tau=60,t=90,a=3. Constant v=0.5, no backreaction.

Native N64 PV2,dt=0.01; frames every 20 steps (dtau=0.2), history
every 10. Full N32 run is a spatial diagnostic only. Refine the same N64
map at dt=0.005 to tau=50; require final theta RMS difference <0.002.
Independent N8 DOP853 of the same cosine semi-discrete PDE to tau=60,
rtol=2e-10,atol=2e-12,max_step=0.02; require max final theta error <0.003,
max pi error <0.03. Require exact map transfer and each frame's mass/scale
agreement. The lower H/m allows more small-scale field structure: N64 is
a representative animation resolution, not a certified continuum result.
F analytic scale identity <1e-12. Do not alter tolerances after execution.
''')
    records=[];states={};gates=[]
    def gate(name,value,limit):gates.append(dict(name=name,value=float(value),limit=limit,passed=bool(np.isfinite(value) and value<=limit)))
    def execute(n,dt,duration,name,frames=False):
        # Scale the already sampled normalized map exactly. New N8 has its own
        # exact periodic projection from the same rescaled catalogue.
        mapfile=root/f'N{n}.afopt'
        if not mapfile.exists():
            if n in (32,64):
                _,previous=read_state(ROOT/f'outputs/stochastic/field-reference/N{n}_dt0.005/initial.bin')
                arrival=previous[...,2]/h
            else:arrival=catalogue.lattice_arrival_times(n)
            write_arrival_map(mapfile,arrival,box)
        params=dict(N=n,lSide=box,dt=dt,tMax=duration,evolver='PV2',expansion='true',fixedBackground='true',H0=h,
                    omegaEoS=1/3,tau_origin=origin,omega_star=1.,f_axion=1.,theta_initial=1.,theta_prime_initial=0.,
                    mass_false=0.,mass_true=1.,potential='cosine',wall_width=width,profile_wall_speed=.5,
                    arrival_map=str(mapfile),result_dir=str(root/name),snapshot_format='binary',output_every=round(.1/dt),profile='true')
        if frames:params['snapshot_every']=20
        (root/('input.json' if frames else f'{name}.input.json')).write_text(json.dumps(params,indent=2)+'\n')
        print(f'EVOLVE {name}: N={n}, {round(duration/dt)} steps',flush=True)
        begin=time.perf_counter()
        with (root/f'{name}.log').open('w') as log:
            subprocess.run([str(ROOT/'build/cosmolattice/axion_fopt_evolve'),*[f'{k}={v}' for k,v in params.items()]],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1200)
        elapsed=time.perf_counter()-begin
        meta,final=read_state(root/name/'final.bin');_,initial=read_state(root/name/'initial.bin')
        np.testing.assert_array_equal(final[...,2],initial[...,2])
        gate(name+'/scale',abs(meta['scale']-h*(origin+duration)),1e-11)
        expected=.5*(1+np.tanh((origin+duration-final[...,2])*.5/width))
        gate(name+'/mass',np.max(np.abs(final[...,3]-expected)),3e-14)
        records.append(dict(name=name,n=n,dt=dt,duration=duration,seconds=elapsed,
                            theta_mean=float(final[...,0].mean()),theta_std=float(final[...,0].std()),
                            timings=json.loads((root/name/'timings.json').read_text())))
        states[name]=final
        print(f'DONE {name}: {elapsed:.2f}s',flush=True)
    execute(64,.01,40.,'evolution',True)
    execute(64,.005,30.,'refined')
    _,same_time=read_state(root/'evolution/frame_0003000.bin')
    gate('N64/timestep_theta_rms',np.sqrt(np.mean((same_time[...,0]-states['refined'][...,0])**2)),.002)
    execute(32,.01,40.,'N32_diagnostic')
    fine=states['evolution'][...,0].reshape(32,2,32,2,32,2).mean(axis=(1,3,5))
    space=float(np.sqrt(np.mean((fine-states['N32_diagnostic'][...,0])**2)))
    execute(8,.01,40.,'N8_reference')
    arrival=states['N8_reference'][...,2];count=8**3;dx=box/8
    def rhs(u,state):
        a=1+h*u;theta=state[:count].reshape(8,8,8)
        lap=sum(np.roll(theta,1,axis)+np.roll(theta,-1,axis)-2*theta for axis in range(3))/dx**2
        m2=.5*(1+np.tanh((origin+u-arrival)*.5/width))
        return np.r_[state[count:]/a**2,(a**2*lap-a**4*m2*np.sin(theta)).ravel()]
    sol=solve_ivp(rhs,(0,40),np.r_[np.ones(count),np.zeros(count)],method='DOP853',rtol=2e-10,atol=2e-12,max_step=.02)
    if not sol.success:raise RuntimeError(sol.message)
    final=states['N8_reference']
    gate('N8/reference_theta',np.max(np.abs(final[...,0].ravel()-sol.y[:count,-1])),.003)
    gate('N8/reference_pi',np.max(np.abs(final[...,1].ravel()-sol.y[count:,-1])),.03)
    report=dict(status='passed' if all(g['passed'] for g in gates) else 'failed',parameters=cfg,runs=records,gates=gates,
                spatial_theta_rms_diagnostic=space,catalogue='Validated H=1 seed-0 catalogue rescaled exactly in coordinates/time; field independently evolved.',
                limits='Representative resolution only; spatial convergence not established. Frozen birth catalogue ends at tau=50.',
                executable_sha256=hashlib.sha256((ROOT/'build/cosmolattice/axion_fopt_evolve').read_bytes()).hexdigest())
    (root/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    archive=root/'sources'
    for file in (Path(__file__),Path(__file__).with_name('make_gifs.py'),ROOT/'adapters/cosmolattice/evolve.cpp'):
        dest=archive/file.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(file,dest)
    print(json.dumps(report['gates'],indent=2),flush=True)
    if report['status']!='passed':raise AssertionError('Validation failed; inspect validation.json')


if __name__=='__main__':main()
