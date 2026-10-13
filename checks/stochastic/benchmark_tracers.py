"""Isolated warm-JIT tracer scaling baseline; see SCALING_EXPERIMENT.md."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import numpy as np
from run_survival import ROOT, curve
from analytic import Case, survival
from tracer_reference import DEFAULT_SOURCE, load_paper, simulate


def instrument(paper):
    measurements=dict(initialization_seconds=0.,engulfment_seconds=0.,growth_seconds=0.,steps=[])
    base_tracer=paper['NucleationTracerManager']; base_bubble=paper['OptimizedBubbleManager']
    class Tracer(base_tracer):
        def __init__(self,*args,**kwargs):
            begin=time.perf_counter(); super().__init__(*args,**kwargs)
            measurements['initialization_seconds']+=time.perf_counter()-begin
            measurements['initial_tracer_array_bytes']=self.tracers.nbytes+self.active_indices.nbytes
        def update_state(self,centres,radii):
            active=self.num_active_tracers
            begin=time.perf_counter(); super().update_state(centres,radii)
            elapsed=time.perf_counter()-begin
            measurements['engulfment_seconds']+=elapsed
            measurements['steps'].append(dict(active_before=active,bubbles=len(centres),
                                                survivors=self.num_active_tracers,seconds=elapsed,
                                                pair_upper_bound=active*len(centres)))
    class Bubble(base_bubble):
        def __init__(self,*args,**kwargs):
            begin=time.perf_counter(); super().__init__(*args,**kwargs)
            measurements['initialization_seconds']+=time.perf_counter()-begin
        def evolve_bubbles(self,tau):
            begin=time.perf_counter(); super().evolve_bubbles(tau)
            measurements['growth_seconds']+=time.perf_counter()-begin
    return dict(paper,NucleationTracerManager=Tracer,OptimizedBubbleManager=Bubble),measurements


def worker(args):
    paper,provenance=load_paper(args.paper_source)
    import numba
    numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
    case=Case(f'L{args.box:g}',args.box,.5,2.5)
    times=np.linspace(0,2.5,26)
    # Same compiled signatures as production; first pass covers JIT compilation.
    warm_case=Case('warmup',3.,.5,2.5)
    ordinary,internal=simulate(paper,warm_case,92843,512,.05,times)
    profiled,_=instrument(paper)
    check,check_internal=simulate(profiled,warm_case,92843,512,.05,times)
    if not (np.array_equal(ordinary.times,check.times) and np.array_equal(ordinary.positions,check.positions)
            and np.array_equal(internal,check_internal)):
        raise AssertionError('Profiling changed reference results')
    baseline_peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    runs=[]
    for seed in (0,1,2):
        profiled,m=instrument(paper)
        begin=time.perf_counter()
        catalogue,internal=simulate(profiled,case,seed,args.ntracers,.05,times)
        total=time.perf_counter()-begin
        measured=sum(m[key] for key in ('initialization_seconds','engulfment_seconds','growth_seconds'))
        if measured>total or not np.isfinite(total): raise AssertionError('Invalid phase timings')
        m['other_generation_seconds']=total-measured
        probes=np.random.default_rng(np.random.SeedSequence([92844,seed])).uniform(0,case.box,(2048,3))
        begin=time.perf_counter(); arrival=catalogue.arrival_times(probes)
        projection=time.perf_counter()-begin
        observed=curve(arrival,times)
        for values in (observed,internal):
            if np.any(np.diff(values)>0) or values.min()<0 or values.max()>1: raise AssertionError('Invalid survival curve')
        runs.append(dict(seed=seed,generation_seconds=total,probe_arrival_seconds=projection,
                         bubbles=len(catalogue.times),geometric_F=observed.tolist(),internal_F=internal.tolist(),**m))
        np.savez_compressed(args.output/f'seed{seed}.npz',times=catalogue.times,positions=catalogue.positions,
                            geometric_F=observed,internal_F=internal)
    rss_unit=1 if sys.platform=='darwin' else 1024
    result=dict(case=asdict(case),ntracers=args.ntracers,runs=runs,tau=times.tolist(),analytic_F=survival(case,times).tolist(),
                paper_source=provenance,numba_threads=numba.get_num_threads(),numpy=np.__version__,
                warmup_equivalence=True,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*rss_unit,
                warmup_peak_rss_bytes=baseline_peak*rss_unit)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/stochastic/tracer-scaling-reference')
    parser.add_argument('--paper-source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--worker',action='store_true')
    parser.add_argument('--density-match-only',action='store_true',help='Only L=6,N=800000, same three seeds')
    parser.add_argument('--box',type=float)
    parser.add_argument('--ntracers',type=int)
    args=parser.parse_args()
    if args.worker:
        worker(args); return
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copy2(args.paper_source,root/'paper_axion_sim.py.txt')
    results=[]
    for box in ((6.,) if args.density_match_only else (3.,6.)):
        for count in ((800000,) if args.density_match_only else (2048,8192,32768,100000)):
            dest=root/f'L{box:g}_N{count}';dest.mkdir()
            print(f'START L={box:g}, tracers={count}',flush=True)
            subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker','--box',str(box),
                            '--ntracers',str(count),'--output',str(dest),'--paper-source',str(args.paper_source)],check=True)
            result=json.loads((dest/'summary.json').read_text());results.append(result)
            secs=[r['generation_seconds'] for r in result['runs']]
            print(f'DONE L={box:g}, tracers={count}: median {np.median(secs):.3f}s, range {min(secs):.3f}..{max(secs):.3f}s',flush=True)
            (root/'partial.json').write_text(json.dumps(results,indent=2)+'\n')
    files=[Path(__file__),Path(__file__).with_name('tracer_reference.py'),Path(__file__).with_name('analytic.py'),
           Path(__file__).with_name('SCALING_EXPERIMENT.md'),ROOT/'src/axion_fopt_next/core.py']
    report=dict(status='passed',configurations=results,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                rng='tracers: legacy MT19937 seeds 0,1,2; probes: default_rng PCG64 SeedSequence([92844,seed])',
                limits='Three realizations/configuration; cost baseline, not accuracy certification. Process RSS includes JIT/runtime. No field evolution.')
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    archive=root/'sources'
    for name in report['source_sha256']:
        target=archive/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    plot(root,results)


def plot(root,results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    for box,color in ((3.,'#4477aa'),(6.,'#cc6677')):
        selected=[r for r in results if r['case']['box']==box]
        if not selected:
            continue
        n=[r['ntracers'] for r in selected]
        for ax,key in zip(axes,('generation_seconds','engulfment_seconds')):
            vals=np.array([[run[key] for run in r['runs']] for r in selected])
            median=np.median(vals,axis=1)
            ax.errorbar(n,median,yerr=[median-vals.min(axis=1),vals.max(axis=1)-median],color=color,marker='o',capsize=3,label=f'H* L={box:g}')
            ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Number of tracers');ax.set_ylabel('Seconds per realization');ax.grid(alpha=.2)
    axes[0].set_title('Complete catalogue generation');axes[1].set_title('Tracer engulfment (nested component)')
    axes[0].legend()
    fig.suptitle('Paper tracer scaling · warm JIT · 8 Numba threads\nMedian and range of 3 seeds; 50 geometry steps, H*=beta=1')
    for ext in ('png','pdf'):fig.savefig(root/f'tracer_scaling.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
