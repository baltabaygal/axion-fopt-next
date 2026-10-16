"""Modest-grid core-growth workflow test, not a post-collapse abundance experiment."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from model import (ROOT,EXECUTABLE,config,run,read_state,gaussian,cool,defects,energies,radial_coupling)

CONTRACT=dict(source='arXiv:2410.18064v2 IV: s=-1 core growth, then s=1 physical coupling',
    lambda_physical=.04,tau_start=8.,tau_cg=12.,tau_end=14.,L=12.,ell=2.,seed=0,
    chi_star=0.,diffusion_duration=2.,diffusion_dt=.01,
    preparation='implicit lattice diffusion at lambda_diff=lambda_phys*(tau_cg/tau_start)^4; then compiled radiation/PV2 core growth',
    matching='N32 initial state is the block average of frozen cooled N64 data',
    limitations='reduced normalization/times; no mass turn-on; no scaling or yield comparison',
    gates=dict(coupling_error=2e-13,width_growth_error=2e-13,minimum_core_cells=4.,
               energy_snapshot_error=2e-12,relative_energy_balance=.003,
               spatial_field_rms=.05,temporal_field_rms=.001))


def measure(directory,cfg):
    h=np.genfromtxt(directory/'native/history.csv',delimiter=',',names=True)
    expected=np.array([radial_coupling(cfg,u) for u in h['u']])
    dx=cfg['lSide']/cfg['N'];width=1/(h['a']*np.sqrt(2*h['lambda']))
    before=h['tau']<=cfg['core_growth_until']
    width_per_tau=width[before]/h['tau'][before]
    width_error=float(np.max(abs(width_per_tau/width_per_tau[0]-1)))
    meta,final=read_state(directory/'native/final.bin')
    e=energies(final,cfg['lSide'],meta['a'],cfg['lambda'],meta['chi'])
    energy_error=max(abs(e[k]-h[k][-1]) for k in e)
    balance=abs(h['total'][-1]-h['total'][0]-np.trapezoid(h['source']+h['dilution'],h['u']))/max(h['total'][0],h['total'][-1])
    checks=dict(coupling=float(np.max(abs(h['lambda']-expected)))<CONTRACT['gates']['coupling_error'],
                width_linear=width_error<CONTRACT['gates']['width_growth_error'],
                core_resolution=float(np.min(width/dx))>=CONTRACT['gates']['minimum_core_cells'],
                energy=energy_error<CONTRACT['gates']['energy_snapshot_error'],
                balance=balance<CONTRACT['gates']['relative_energy_balance'])
    snapshots=[]
    for path in sorted((directory/'native').glob('frame_*.bin')):
        meta,state=read_state(path);d,_=defects(state,cfg['lSide'],meta['tau'])
        snapshots.append(dict(**meta,**d,lambda_actual=radial_coupling(cfg,meta['u'])))
    report=dict(parameters=cfg,checks={k:bool(v) for k,v in checks.items()},passed=all(checks.values()),
                width_linear_error=width_error,minimum_core_cells=float(np.min(width/dx)),
                balance=float(balance),snapshot_energy_error=float(energy_error),snapshots=snapshots)
    (directory/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return report,final


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=EXECUTABLE)
    p.add_argument('--box',type=float,choices=(12.,24.),default=12.)
    args=p.parse_args()
    root=args.output;root.mkdir(parents=True,exist_ok=False)
    contract=dict(CONTRACT,L=args.box,
        matching=CONTRACT['matching'] if args.box==12 else 'Frozen N64 initial state for timestep comparison; spatial validation is the separate L12 suite',
        require_survival_at_cg=args.box==24)
    (root/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    raw=gaussian(64,args.box,2.,0)
    diffusion_lambda=.04*(12./8.)**4
    initial,cooling=cool(raw,args.box,diffusion_lambda,duration=2.,dt=.01)
    np.save(root/'initial_raw.npy',raw);np.save(root/'initial_cooled.npy',initial)
    np.savetxt(root/'cooling.csv',cooling,delimiter=',',header='diffusion_time,total',comments='')
    coarse=initial.reshape(32,2,32,2,32,2,4).mean(axis=(1,3,5))
    reports={};finals=[]
    plans=[('N32',32,.01,coarse),('N64',64,.01,initial),('N32_fine_dt',32,.005,coarse)] if args.box==12 else [('N64',64,.01,initial),('N64_fine_dt',64,.005,initial)]
    for label,n,dt,state in plans:
        cfg=config(n,args.box,dt,6.,lam=.04,origin=8.,star=20.,chi_star=0.,every=1,frames=round(.5/dt))
        cfg['core_growth_until']=12.
        cfg=run(state,cfg,root/label,args.executable)
        report,final=measure(root/label,cfg);reports[label]=report;finals.append(final)
        print(label,report['passed'],'core minimum',report['minimum_core_cells'],flush=True)
    spatial=None
    if args.box==12:
        blocked=finals[1].reshape(32,2,32,2,32,2,4).mean(axis=(1,3,5))
        spatial=float(np.sqrt(np.mean((finals[0][...,:2]-blocked[...,:2])**2)))
        fine=finals[2]
    else:
        fine=finals[1]
    temporal=float(np.sqrt(np.mean((finals[0][...,:2]-fine[...,:2])**2)))
    checks=dict(temporal=temporal<CONTRACT['gates']['temporal_field_rms'],
                monotone_cooling=bool(np.max(np.diff(cooling[:,1]))<1e-12))
    if spatial is not None:checks['spatial']=spatial<CONTRACT['gates']['spatial_field_rms']
    if args.box==24:
        cg_state=next(r for r in reports['N64']['snapshots'] if r['tau']==12.)
        checks['network_survives_at_cg']=cg_state['plaquettes']>0
    passed=all(r['passed'] for r in reports.values()) and all(checks.values())
    result=dict(contract=contract,runs=reports,comparisons=dict(checks=checks,spatial_rms=spatial,temporal_rms=temporal),
                passed=passed,executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'checks/pq_strings/pq_model.h',ROOT/'checks/pq_strings/evolve.cpp',ROOT/'checks/pq_strings/model.py')})
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Preparation suite passed:',passed);raise SystemExit(0 if passed else 1)
if __name__=='__main__':main()
