"""Evolved-wave validation of the current estimator through rapid mass growth.

A defect-free tiny-amplitude Fourier mode is a cheap estimator test even when
one cannot resolve a full string network. It is not a network-abundance result.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import EXECUTABLE,config,run,read_state,coefficients
from spectra import currents,mode_spectrum,normalized_yield

GATES=dict(max_relative_yield_error=.005,max_angle_error_over_initial=.001,
           timestep_yield_change=.001,late_plateau_peak_to_peak=.05,
           minimum_late_mass_phase_advance=4*np.pi)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=EXECUTABLE);args=p.parse_args()
    root=args.output;root.mkdir(parents=True,exist_ok=False)
    (root/'contract.json').write_text(json.dumps(dict(gates=GATES,
        initial='N8,L8,theta=1e-5*cos(2pi*x/L); |Phi|=1; WKB initial angular velocity',
        physics='lambda=.08,tau_origin=20,tau_star=27,n_tau=8; evolve to tau=81',
        limitation='single resolved linear mode, no defects; not a network reproduction'),indent=2)+'\n')
    n=8;box=8.;amplitude=1e-5;k=2*np.sin(np.pi/n)/(box/n)
    x=(np.arange(n)+.5)*box/n;cosine=np.broadcast_to(np.cos(2*np.pi*x/box)[:,None,None],(n,n,n))
    cfg=config(n,box,.005,61.,every=100,frames=100)
    _,a0,chi0=coefficients(cfg,0.)
    omega0=np.sqrt(k*k+a0*a0*chi0)
    omega_prime0=5*a0*a0*chi0/(cfg['tau_origin']*omega0)
    angular_velocity=-(1/cfg['tau_origin']+omega_prime0/(2*omega0))*amplitude
    theta=amplitude*cosine;phi=np.exp(1j*theta);momentum=a0*a0*1j*phi*angular_velocity*cosine
    state=np.stack((phi.real,phi.imag,momentum.real,momentum.imag),axis=-1)
    def rhs(u,y):
        _,a,chi=coefficients(cfg,u)
        return [y[1]/a**2,-a*a*k*k*y[0]-a**4*chi*y[0]]
    exact=solve_ivp(rhs,(0.,61.),[amplitude,a0*a0*angular_velocity],method='DOP853',
                    rtol=2e-11,atol=2e-16,max_step=.01,dense_output=True)
    if not exact.success:raise RuntimeError(exact.message)
    checks={};details={};histories=[]
    for dt in (.005,.0025):
        cfg['dt']=dt;cfg['snapshot_every']=round(.5/dt);cfg['output_every']=round(.5/dt)
        directory=root/f'dt{dt}';actual=run(state,cfg,directory,args.executable)
        rows=[];yield_errors=[];angle_errors=[]
        for path in sorted((directory/'native').glob('frame_*.bin')):
            meta,evolved=read_state(path);_,a,chi=coefficients(cfg,meta['u'])
            j0,ji,_=currents(evolved,box,a)
            modes=mode_spectrum(j0,ji,box,a,np.sqrt(chi),dispersion='lattice')
            th,p=exact.sol(meta['u']);omega=np.sqrt(k*k+a*a*chi)
            expected=(p*p/a**4/omega+omega*th*th)/(4*a)
            yield_errors.append(abs(modes['summary']['number_density']/expected-1))
            angle_errors.append(float(np.max(abs(np.arctan2(evolved[...,1],evolved[...,0])-th*cosine)))/amplitude)
            nu=normalized_yield(modes['summary']['number_density'],a,27/20,20/27**2)
            rows.append(dict(u=meta['u'],tau_ratio=meta['tau']/27,number=modes['summary']['number_density'],
                             expected_number=float(expected),nu=nu,omega=float(omega)))
        late=[r for r in rows if r['tau_ratio']>=2.7]
        spread=(max(r['nu'] for r in late)-min(r['nu'] for r in late))/np.mean([r['nu'] for r in late])
        phase=float(np.trapezoid([r['omega'] for r in late],[r['u'] for r in late]))
        checks[f'{dt}_yield']=max(yield_errors)<GATES['max_relative_yield_error']
        checks[f'{dt}_field']=max(angle_errors)<GATES['max_angle_error_over_initial']
        checks[f'{dt}_plateau']=spread<GATES['late_plateau_peak_to_peak']
        checks[f'{dt}_phase']=phase>=GATES['minimum_late_mass_phase_advance']
        details[str(dt)]=dict(max_relative_yield_error=max(yield_errors),max_angle_error=max(angle_errors),
                             late_peak_to_peak=float(spread),late_phase_advance=phase,history=rows)
        histories.append(np.array([r['nu'] for r in rows]))
    change=float(np.max(abs(histories[0]/histories[1]-1)))
    checks['timestep']=change<GATES['timestep_yield_change']
    report=dict(gates=GATES,checks={k:bool(v) for k,v in checks.items()},details=details,
                timestep_relative_change=change,passed=all(checks.values()),
                interpretation='validated estimator and evolved isolated-mode adiabatic yield, not a cosmological network yield')
    (root/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='details'},indent=2))
    raise SystemExit(0 if report['passed'] else 1)
if __name__=='__main__':main()
