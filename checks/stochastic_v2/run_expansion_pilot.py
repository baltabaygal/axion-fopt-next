"""Controlled expanding single bubble, independent diagnostics and refinements."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from model import seed,read_state
from expanding_model import invoke,physical_terms,balances,scale,cosmic
from run_bubble import contour
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--validation',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();build=verify_build(args.executable);v=json.loads(args.validation.read_text())
    if not v['passed'] or v['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Matching reference validation required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);shutil.copy2(Path(__file__).with_name('EXPANSION_STOCHASTIC_CONTRACT.md'),out/'contract.md')
    sources=out/'sources';sources.mkdir()
    for source in Path(__file__).parent.iterdir():
        if source.suffix in ('.py','.cpp','.h','.md'):shutil.copy2(source,sources/source.name)
    h0=.025;box=24.;end=16.;results={};timings=[]
    for name,n,dt in [('main',96,.02),('coarse',64,.02),('time-half',96,.01)]:
        dest=out/name;seconds=invoke(args.executable,dest,seed(n,box),box,dt,end,h0)
        h=np.genfromtxt(dest/'native/history.csv',delimiter=',',names=True);_,final=read_state(dest/'native/final.bin');_,initial=read_state(dest/'native/initial.bin')
        a=scale(end,h0);energies=physical_terms(final,box/n,a);error=max(abs(float(value.sum()*(box/n)**3)-h[key][-1]) for key,value in energies.items())/h['total'][0]
        total,axion=balances(h);wall=[];stop=None
        for f in sorted((dest/'native').glob('frame_*.bin')):
            meta,plane=read_state(f);c=contour(plane[:,n//2-1:n//2+1,0].mean(axis=1),box,(box/n)**2/2)
            if c['width'] is None or c['r10']>=box/2:stop=meta['t'];break
            wall.append(dict(u=meta['t'],a=float(scale(meta['t'],h0)),comoving_width=c['width'],physical_width=float(scale(meta['t'],h0)*c['width']),cells=c['width']/(box/n)))
        checks=dict(initial_exact=bool(np.array_equal(initial,seed(n,box))),finite=bool(np.isfinite(final).all()),snapshot_energy=error<1e-11,
            total_balance=np.max(abs(total))/h['total'][0]<.005,axion_balance=np.max(abs(axion))/h['total'][0]<.005,
            expansion_transfer=np.max(abs(h['a']-scale(h['u'],h0)))<1e-12 and np.max(abs(h['t_cosmic']-cosmic(h['u'],h0)))<1e-12)
        results[name]=dict(checks={k:bool(value) for k,value in checks.items()},n=n,dt=dt,box=box,H0=h0,u_final=end,cosmic_elapsed_final=float(cosmic(end,h0)),
            a_final=float(a),axion_energy_final=float(h['E_A'][-1]),phase_fraction_final=float(1-h['F_half'][-1]),
            total_balance=float(np.max(abs(total))/h['total'][0]),axion_balance=float(np.max(abs(axion))/h['total'][0]),snapshot_energy_error=float(error),
            wall_resolution=wall,radial_diagnostic_stop=stop)
        timings.append(dict(case=name,seconds=seconds));(out/'timings.json').write_text(json.dumps(timings,indent=2)+'\n')
        print(json.dumps(dict(case=name,seconds=seconds,checks=results[name]['checks'])),flush=True)
    main=results['main'];spatial=abs(results['coarse']['axion_energy_final']/main['axion_energy_final']-1);temporal=abs(results['time-half']['axion_energy_final']/main['axion_energy_final']-1)
    checks=dict(local=all(all(r['checks'].values()) for r in results.values()),spatial_energy=spatial<.05,temporal_energy=temporal<.01,
        resolved_preinteraction_wall=bool(main['wall_resolution']) and min(r['cells'] for r in main['wall_resolution'])>=5)
    report=dict(passed=bool(all(checks.values())),checks=checks,results=results,spatial_energy_difference=spatial,temporal_energy_difference=temporal,
        build=build,validation_sha256=sha(args.validation),scope='Prescribed RD, polynomial coupled fields; single seeded periodic bubble; no Friedmann, reheating or abundance claim')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('build','results')},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
