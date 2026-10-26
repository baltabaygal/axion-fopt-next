"""Extend a frozen stochastic realization and distinguish conversion/relaxation."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from scipy.integrate import quad
from model import read_state
from expanding_model import invoke,scale,cosmic,balances
from stochastic_reference import insertion,energy_totals,rate
from validate_backend import verify_build
from build_backend import sha

def intervals(h,mask):
    # Use the post-event value where multiple histories have the same time.
    keep=np.r_[np.diff(h['u'])>1e-10,True];u=h['u'][keep];mask=mask[keep]
    spans=[];start=None
    for i,yes in enumerate(mask):
        if yes and start is None:start=float(u[i])
        if start is not None and (not yes or i==len(u)-1):
            end=float(u[i-1] if not yes else u[i]);spans.append(dict(start=start,end=end,duration=end-start));start=None
    return spans

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('executable','validation','candidates','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--case',choices=['main','space','time-half'],default='main')
    p.add_argument('--ordinary-fine',action='store_true',help='Keep the late timestep cap but disable certified rejection')
    p.add_argument('--end',type=float,help='Continue beyond catalogue only with global empty-eligibility certificates, maximum 80')
    p.add_argument('--stop-sustained',type=float,default=0,help='Stop after this all-cell interval, 0 disables')
    p.add_argument('--previous',type=Path,default=Path('outputs/stochastic_v2/expanding-stochastic-01'))
    args=p.parse_args();build=verify_build(args.executable);validation=json.loads(args.validation.read_text())
    audit=json.loads((args.candidates/'report.json').read_text());end=args.end or audit['end']
    if end<audit['end'] or end>80:raise ValueError('Require catalogue endpoint <= end <= 80')
    planned_end=end
    if not validation['passed'] or validation['build']['executable_sha256']!=build['executable_sha256'] or not audit['passed']:
        raise RuntimeError('Matching validation and passing candidate audit required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);sources=out/'sources';sources.mkdir()
    for source in Path(__file__).parent.iterdir():
        if source.suffix in ('.py','.cpp','.h','.md'):shutil.copy2(source,sources/source.name)
    shutil.copy2(args.candidates/'candidates.txt',out/'candidates.txt');shutil.copy2(args.candidates/'report.json',out/'candidate-audit.json')
    n=128 if args.case=='space' else 96;dt=.01 if args.case=='time-half' else .02
    box=32.;h0=.025;y=np.zeros((n,n,n,4));y[...,1]=.2
    extra=dict(candidate_file=str(out/'candidates.txt'),seed_radius=8,seed_width=2,core_tolerance=.05)
    if args.ordinary_fine:extra['certify_rejections']='false'
    if end>audit['end']:extra['certify_empty_after']=audit['end']
    if args.stop_sustained:extra['stop_sustained']=args.stop_sustained
    dest=out/args.case;seconds=invoke(args.executable,dest,y,box,dt,end,h0,extra=extra,timeout=1800)
    native=dest/'native';h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
    ledger=np.atleast_1d(np.genfromtxt(native/'births.csv',delimiter=',',names=True));final_meta,final=read_state(native/'final.bin');_,initial=read_state(native/'initial.bin')
    end=final_meta['t']
    normal=max(abs(h['total'][0]),abs(h['total'][-1]),float(sum(abs(ledger['delta_total']))))
    et=energy_totals(final,box,scale(end,h0));endpoint=float(max(abs(et-np.array([h['total'][-1],h['E_S'][-1],h['E_A'][-1]])))/normal)
    total,axion=balances(h);maps=[];jumps=[];unchanged=[]
    for i,row in enumerate(ledger):
        if not row['accepted']:continue
        meta,pre=read_state(native/f'birth_pre_{i+1}.bin');_,post=read_state(native/f'birth_post_{i+1}.bin')
        expected,decision=insertion(pre,box,meta['t'],h0,[row['x'],row['y'],row['z']]);assert decision['accepted']
        maps.append(float(max(abs(post-expected).ravel())));unchanged.append(np.array_equal(post[...,1::2],pre[...,1::2]))
        jump=energy_totals(post,box,scale(meta['t'],h0))-energy_totals(pre,box,scale(meta['t'],h0))
        jumps.append(float(max(abs(jump-np.array([row['delta_total'],row['delta_S'],row['delta_A']])))/max(1,float(max(abs(jump))))))
    overlap=None
    if args.case!='space':
        old=np.genfromtxt(args.previous/args.case/'native/history.csv',delimiter=',',names=True)
        prefix=h[h['u']<=old['u'][-1]+1e-10]
        if len(prefix)!=len(old):raise RuntimeError('Prefix output scheduling differs')
        overlap=float(max(max(abs(prefix[k]-old[k])) for k in old.dtype.names))
    s=final[...,0];phase_error=abs(float(np.mean(s<.5))-h['F_half'][-1])
    stat_error=max(abs(float(s.min())-h['S_min'][-1]),abs(float(s.max())-h['S_max'][-1]),
        abs(float(s.mean())-h['S_mean'][-1]),abs(float(s.std())-h['S_std'][-1]))
    converted=intervals(h,h['S_min']>=.5);relaxed=intervals(h,(h['S_min']>.9)&(h['S_max']<1.1))
    longest=max([v['duration'] for v in converted],default=0);accepted=ledger['accepted'].astype(bool)
    checks=dict(initial_exact=np.array_equal(y,initial),finite=np.isfinite(final).all(),snapshot_energy=endpoint<1e-11,
        total_balance=max(abs(total))/normal<.005,axion_balance=max(abs(axion))/normal<.005,
        maps=bool(maps) and max(maps)<2e-14,jumps=bool(jumps) and max(jumps)<1e-10,axion_unchanged=all(unchanged),
        endpoint_phase=phase_error==0 and stat_error<1e-12,prefix_overlap=overlap is None or overlap<1e-10)
    empty=None
    if end>audit['end']:
        proof=np.atleast_1d(np.genfromtxt(native/'empty_eligibility.csv',delimiter=',',names=True))
        checks['empty_eligibility']=bool(len(proof)>0 and abs(proof['u_start'][0]-audit['end'])<1e-10 and
            abs(proof['u_end'][-1]-end)<1e-10 and np.max(abs(proof['u_start'][1:]-proof['u_end'][:-1]),initial=0)<1e-10 and
            np.all(proof['u_end']>proof['u_start']) and np.all(proof['min_witness']>.05) and np.all(proof['envelope']<1))
        empty=dict(start=audit['end'],end=end,steps=len(proof),minimum_witness=float(proof['min_witness'].min()),
            expected_unmaterialized_attempts=float(box**3*quad(rate,audit['end'],end)[0]),
            interpretation='No realized count/positions drawn; every possible centre analytically rejected throughout every covered step')
    report=dict(passed=bool(all(checks.values())),checks={k:bool(v) for k,v in checks.items()},case=args.case,n=n,dt=dt,box=box,H0=h0,
        u_final=end,cosmic_elapsed_final=float(cosmic(end,h0)),a_final=float(scale(end,h0)),seconds=seconds,
        candidate_count=len(ledger),accepted_count=int(accepted.sum()),accepted_indices=np.flatnonzero(accepted).tolist(),
        certified_rejections=int(ledger['certified_rejection'].sum()) if 'certified_rejection' in ledger.dtype.names else 0,
        ordinary_fine=args.ordinary_fine,
        requested_end=planned_end,stop_sustained=args.stop_sustained,
        empty_eligibility=empty,
        maximum_late_dt=dt*.25 if 'certified_rejection' in ledger.dtype.names else dt,
        accepted_times=ledger['u'][accepted].tolist(),all_cell_intervals=converted,relaxed_intervals=relaxed,
        sustained_conversion=longest>=2-1e-10,longest_converted_interval=longest,
        phase_fraction_final=float(1-h['F_half'][-1]),near_true_fraction_final=float(h['near_true_fraction'][-1]),
        S_min_final=float(s.min()),S_max_final=float(s.max()),S_std_final=float(s.std()),
        axion_energy_final=float(h['E_A'][-1]),total_balance=float(max(abs(total))/normal),axion_balance=float(max(abs(axion))/normal),
        energy_normalization=float(normal),snapshot_energy_error=endpoint,endpoint_stat_error=float(stat_error),
        prefix_overlap_error=overlap,max_map_error=max(maps),max_jump_error=max(jumps),
        seed_width_cells=(2*np.arctanh(.8)*2/(scale(ledger['u'][accepted],h0)*box/n)).tolist(),
        build=build,validation_sha256=sha(args.validation),candidate_sha256=sha(out/'candidates.txt'),
        script_sha256=sha(__file__),scope='One conditional Poisson extension, phase thresholds distinct from vacuum relaxation')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='build'},indent=2),flush=True)
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
