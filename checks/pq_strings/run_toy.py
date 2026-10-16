"""Finite periodic-box capped-mass toy: explicitly excludes paper reproduction."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from model import config,run,deterministic,gaussian,cool,read_state,defects,radial_coupling,ROOT,EXECUTABLE

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--end',type=float,default=80.);p.add_argument('--n',type=int,default=64)
    p.add_argument('--dt',type=float,default=.02);p.add_argument('--seed',type=int)
    p.add_argument('--control',action='store_true');p.add_argument('--frozen',type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    origin=40.;cg=44.;ramp=6.;cap=.06;lam=.08;box=24.
    def chi(t):
        s=np.clip((t-cg)/ramp,0.,1.);return cap**2*s**3*(10+s*(-15+6*s))
    star=brentq(lambda t:np.sqrt(chi(t))-origin/t**2,cg,cg+ramp)
    cfg=config(a.n,box,a.dt,a.end-origin,lam=lam,origin=origin,star=star,every=round(.2/a.dt),frames=round(2/a.dt))
    cfg.update(core_growth_until=cg,mass_schedule='smooth_cap',mass_onset=cg,mass_ramp_duration=ramp,mass_cap=0. if a.control else cap)
    contract=dict(interpretation='Modified finite-box toy; no cosmological abundance or paper reproduction',
      parameters=cfg,initial='symmetric four straight strings' if a.seed is None else f'Gaussian ell=2.5 seed={a.seed}',
      gates=dict(core_cells_while_defects=4,wall_cells=4,max_ma_ms=.16,late_theta_max=.3,late_radial_rms=.1,
      late_window_duration=60,late_phase_advance=4*np.pi,late_yield_peak_to_peak=.05,dt_yield_change=.01,spatial_yield_change=.1),
      initial_freeze='Exact N64 state; N128 Fourier interpolation, matched physical parameters',
      normalization='nu_toy=(a/a_star)^3 n_a/H_star, f_a=1, star solves m_a=H for massive run; control uses same star')
    (a.output/'contract.json').write_text(json.dumps(contract,indent=2))
    if a.frozen:
        initial=np.load(a.frozen)
        if len(initial)!=a.n:
            from scipy.signal import resample
            for axis in range(3):initial=resample(initial,a.n,axis=axis)
    else:
        ld=lam*(cg/origin)**4
        initial=deterministic(a.n,box,core=1/np.sqrt(2*ld)) if a.seed is None else gaussian(a.n,box,2.5,a.seed)
        initial,_=cool(initial,box,ld,duration=2,dt=.01)
    np.save(a.output/'frozen.npy',initial)
    provenance=dict(executable_sha256=hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest(),sources={})
    for name in ('pq_model.h','evolve.cpp','model.py','run_toy.py','spectra.py'):
        provenance['sources'][name]=hashlib.sha256((ROOT/'checks/pq_strings'/name).read_bytes()).hexdigest()
    (a.output/'provenance.json').write_text(json.dumps(provenance,indent=2))
    cfg=run(initial,cfg,a.output/'run')
    rows=[]
    for file in sorted((a.output/'run/native').glob('frame_*.bin')):
        meta,s=read_state(file);d,_=defects(s,box,meta['tau']);ms=np.sqrt(2*radial_coupling(cfg,meta['u']));dx=box/a.n
        rows.append(dict(tau=meta['tau'],plaquettes=d['plaquettes'],walls=d['wall_links'],core_cells=1/(meta['a']*ms*dx),theta_max=float(np.max(abs(np.arctan2(s[...,1],s[...,0]))))))
    (a.output/'topology.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
