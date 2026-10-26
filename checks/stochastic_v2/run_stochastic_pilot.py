"""Frozen Poisson candidates, live finite-core rejection, and controlled RD runs."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from model import read_state
from expanding_model import invoke,scale,cosmic,balances
from stochastic_reference import insertion,energy_totals
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','validation','rate-audit','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();build=verify_build(args.executable)
    validation=json.loads(args.validation.read_text());audit=json.loads((args.rate_audit/'report.json').read_text())
    if not validation['passed'] or validation['build']['executable_sha256']!=build['executable_sha256'] or not audit['passed']:
        raise RuntimeError('Matching insertion validation and passing rate audit required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('candidates.txt','rate.csv'):shutil.copy2(args.rate_audit/name,out/name)
    shutil.copy2(Path(__file__).with_name('EXPANSION_STOCHASTIC_CONTRACT.md'),out/'contract.md')
    sources=out/'sources';sources.mkdir()
    for source in Path(__file__).parent.iterdir():
        if source.suffix in ('.py','.cpp','.h','.md'):shutil.copy2(source,sources/source.name)
    box=32.;end=12.;h0=.025;results={};timings=[]
    candidates=np.atleast_2d(np.loadtxt(out/'candidates.txt'))
    for name,n,dt in [('main',96,.02),('coarse',64,.02),('time-half',96,.01)]:
        initial=np.zeros((n,n,n,4));initial[...,1]=.2;dest=out/name
        seconds=invoke(args.executable,dest,initial,box,dt,end,h0,extra=dict(candidate_file=str(out/'candidates.txt'),seed_radius=8,seed_width=2,core_tolerance=.05))
        native=dest/'native';h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
        ledger=np.atleast_1d(np.genfromtxt(native/'births.csv',delimiter=',',names=True));_,final=read_state(native/'final.bin')
        _,saved_initial=read_state(native/'initial.bin');jump_errors=[];map_errors=[];axion_unchanged=[];core_errors=[]
        for i,row in enumerate(ledger):
            if not row['accepted']:continue
            m,pre=read_state(native/f'birth_pre_{i+1}.bin');_,post=read_state(native/f'birth_post_{i+1}.bin')
            expected,decision=insertion(pre,box,m['t'],h0,[row['x'],row['y'],row['z']])
            map_errors.append(float(np.max(abs(post-expected))))
            core_errors.append(abs(decision['maxcore']-row['max_core_abs_S']))
            axion_unchanged.append(np.array_equal(pre[...,1::2],post[...,1::2]))
            jump=energy_totals(post,box,scale(m['t'],h0))-energy_totals(pre,box,scale(m['t'],h0))
            logged=np.array([row['delta_total'],row['delta_S'],row['delta_A']])
            jump_errors.append(float(max(abs(jump-logged))/max(1,float(max(abs(jump))))))
        energies=energy_totals(final,box,scale(end,h0));native_energies=np.array([h['total'][-1],h['E_S'][-1],h['E_A'][-1]])
        normal=max(abs(h['total'][0]),abs(h['total'][-1]),sum(abs(ledger['delta_total'])))
        endpoint=float(max(abs(energies-native_energies))/normal);total,axion=balances(h)
        candidate_error=float(np.max(abs(np.column_stack([ledger[k] for k in ('u','x','y','z')])-candidates)))
        accepted=ledger['accepted'].astype(bool);birth_scale=scale(ledger['u'][accepted],h0)
        # The analytic seed's 10--90% width is 2*atanh(.8)*w; it is not a late-time wall measurement.
        seed_cells=2*np.arctanh(.8)*2/(birth_scale*(box/n))
        checks=dict(initial_exact=np.array_equal(saved_initial,initial),finite=np.isfinite(final).all(),
            frozen_candidates=candidate_error<1e-11,nontrivial=len(map_errors)>1,snapshot_energy=endpoint<1e-11,
            maps=bool(map_errors) and max(map_errors)<2e-14,jumps=bool(jump_errors) and max(jump_errors)<1e-10,
            axion_unchanged=all(axion_unchanged),core_measure=bool(core_errors) and max(core_errors)<1e-14,
            total_balance=max(abs(total))/normal<.005,axion_balance=max(abs(axion))/normal<.005,
            scale_transfer=max(abs(h['a']-scale(h['u'],h0)))<1e-12,
            resolved_seed=bool(len(seed_cells)) and min(seed_cells)>=5)
        results[name]=dict(checks={k:bool(v) for k,v in checks.items()},n=n,dt=dt,box=box,H0=h0,u_final=end,
            cosmic_elapsed_final=float(cosmic(end,h0)),a_final=float(scale(end,h0)),candidate_count=len(ledger),
            accepted=accepted.tolist(),accepted_count=int(accepted.sum()),axion_energy_final=float(h['E_A'][-1]),
            phase_fraction_final=float(1-h['F_half'][-1]),S_min_final=float(h['S_min'][-1]),S_max_final=float(h['S_max'][-1]),
            total_balance=float(max(abs(total))/normal),axion_balance=float(max(abs(axion))/normal),
            snapshot_energy_error=endpoint,max_map_error=max(map_errors,default=0),max_jump_error=max(jump_errors,default=0),
            insertion_work=ledger['delta_total'].tolist(),axion_insertion_work=ledger['delta_A'].tolist(),
            seed_width_cells=seed_cells.tolist(),energy_normalization=float(normal))
        timings.append(dict(case=name,seconds=seconds));(out/'timings.json').write_text(json.dumps(timings,indent=2)+'\n')
        print(json.dumps(dict(case=name,seconds=seconds,result=results[name])),flush=True)
    main=results['main'];same_space=results['coarse']['accepted']==main['accepted'];same_time=results['time-half']['accepted']==main['accepted']
    spatial=abs(results['coarse']['axion_energy_final']/main['axion_energy_final']-1)
    temporal=abs(results['time-half']['axion_energy_final']/main['axion_energy_final']-1)
    checks=dict(local=all(all(r['checks'].values()) for r in results.values()),
        spatial_energy=(spatial<.1 if same_space else True),temporal_energy=(temporal<.02 if same_time else True))
    report=dict(passed=bool(all(checks.values())),checks=checks,results=results,
        same_spatial_acceptance=same_space,same_temporal_acceptance=same_time,
        paired_spatial_convergence=same_space,paired_temporal_convergence=same_time,
        spatial_energy_difference=spatial,temporal_energy_difference=temporal,build=build,
        validation_sha256=sha(args.validation),rate_audit_sha256=sha(args.rate_audit/'report.json'),
        candidate_sha256=sha(out/'candidates.txt'),rate_sha256=sha(out/'rate.csv'),rate_audit=audit,
        scope='One finite-radius seeded realization; prescribed RD; polynomial reciprocal coupling; no bounce, bath, abundance or continuum claim')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('results','build','rate_audit')},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
