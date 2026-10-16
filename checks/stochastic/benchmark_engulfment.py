"""Paired periodic spatial-index comparison against frozen all-pairs kernel."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import numpy as np
from run_survival import ROOT,curve
from tracer_reference import DEFAULT_SOURCE,load_paper,simulate
from analytic import Case
from benchmark_tracers import instrument
from spatial_engulfment import accelerated_paper,build_index,query_index


def verify_kernel(paper):
    rng=np.random.default_rng(92845); oracle=paper['find_surviving_tracers_simple'];checked=0
    cases=[]
    for box in (1.,3.,6.):
        for bubbles in (0,1,10,200):
            for radius in (0.,.03,.3,1.):
                centres=rng.uniform(-box/2,box/2,(bubbles,3)).astype(np.float32)
                radii=rng.uniform(0,radius*box,bubbles).astype(np.float32)
                tracers=rng.uniform(-box/2,box/2,(1000,3)).astype(np.float32)
                cases.append((box,tracers,centres,radii))
    # Exact tangent, nextafter, periodic face/corner, cell boundaries, +/-L/2.
    for box in (1.,3.,6.):
        b=np.float32(box);r=np.float32(box/4)
        centres=np.array([[0,0,0],[box/2,box/2,box/2]],dtype=np.float32)
        radii=np.array([r,np.float32(box/8)],dtype=np.float32)
        points=[]
        for axis in range(3):
            for value in (r,np.nextafter(r,np.float32(0)),np.nextafter(r,np.float32(np.inf)),b/2,-b/2):
                p=np.zeros(3,dtype=np.float32);p[axis]=value;points.append(p)
        points.extend([[box/2,box/2,box/2],[-box/2,-box/2,-box/2]])
        cases.append((box,np.array(points,dtype=np.float32),centres,radii))
        cases.append((box,np.empty((0,3),dtype=np.float32),centres,radii))
    for box,tracers,centres,radii in cases:
        active=rng.permutation(len(tracers)).astype(np.int64)
        reference=oracle(tracers,active,centres,radii,box,8)
        for ncells in (1,4,16):
            offsets,entries=build_index(centres,radii,box,ncells,8_000_000)
            actual=query_index(tracers,active,centres,radii,box,ncells,offsets,entries)
            np.testing.assert_array_equal(actual,reference);checked+=1
    # Force class fallback and compare a complete realization, not just query.
    case=Case('fallback',3.,.5,2.5);times=np.linspace(0,2.5,26)
    a,fa=simulate(paper,case,92846,512,.05,times)
    b,fb=simulate(accelerated_paper(paper,max_references=0),case,92846,512,.05,times)
    np.testing.assert_array_equal(a.times,b.times);np.testing.assert_array_equal(a.positions,b.positions)
    np.testing.assert_array_equal(fa,fb)
    return dict(identical_input_checks=checked,fallback_equivalent=True)


def worker(args):
    paper,provenance=load_paper(args.paper_source)
    import numba
    numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
    validation=verify_kernel(paper)
    optimized=accelerated_paper(paper)
    warm=Case('warmup',3.,.5,2.5);times=np.linspace(0,2.5,26)
    for engine in (paper,optimized):simulate(engine,warm,92843,512,.05,times)
    case=Case(f'L{args.box:g}',args.box,.5,2.5);pairs=[]
    for seed in (0,1,2):
        engines={'baseline':paper,'indexed':optimized};rows={};catalogues={};fractions={}
        order=('baseline','indexed') if seed%2==0 else ('indexed','baseline')
        for method in order:
            profiled,m=instrument(engines[method]);begin=time.perf_counter()
            catalogue,internal=simulate(profiled,case,seed,args.ntracers,.05,times)
            elapsed=time.perf_counter()-begin
            parts=sum(m[k] for k in ('initialization_seconds','engulfment_seconds','growth_seconds'))
            if not np.isfinite(elapsed) or parts>elapsed:raise AssertionError('Bad timings')
            m['other_generation_seconds']=elapsed-parts
            rows[method]=dict(generation_seconds=elapsed,**m)
            catalogues[method]=catalogue;fractions[method]=internal
        a,b=catalogues['baseline'],catalogues['indexed']
        np.testing.assert_array_equal(a.times,b.times);np.testing.assert_array_equal(a.positions,b.positions)
        np.testing.assert_array_equal(fractions['baseline'],fractions['indexed'])
        np.testing.assert_array_equal([s['survivors'] for s in rows['baseline']['steps']],
                                      [s['survivors'] for s in rows['indexed']['steps']])
        probes=np.random.default_rng(np.random.SeedSequence([92844,seed])).uniform(0,case.box,(2048,3))
        arrival_a=a.arrival_times(probes);arrival_b=b.arrival_times(probes)
        np.testing.assert_array_equal(arrival_a,arrival_b)
        rows['seed']=seed;rows['order']=order;rows['bubbles']=len(a.times);rows['exact_equivalence']=True
        pairs.append(rows)
        np.savez_compressed(args.output/f'seed{seed}.npz',times=a.times,positions=a.positions,
                            internal_F=fractions['baseline'],geometric_F=curve(arrival_a,times))
        print(f'  seed {seed}: {rows["baseline"]["generation_seconds"]:.3f} -> {rows["indexed"]["generation_seconds"]:.3f}s (identical)',flush=True)
    result=dict(box=args.box,ntracers=args.ntracers,pairs=pairs,validation=validation,paper_source=provenance,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
                numba_threads=numba.get_num_threads(),numpy=np.__version__,ncells=16,max_index_references=8_000_000)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/stochastic/engulfment-index-reference')
    parser.add_argument('--paper-source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--box',type=float);parser.add_argument('--ntracers',type=int)
    args=parser.parse_args()
    if args.worker:worker(args);return
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copy2(args.paper_source,root/'paper_axion_sim.py.txt');results=[]
    for box,count in ((3.,100000),(6.,100000),(6.,800000)):
        dest=root/f'L{box:g}_N{count}';dest.mkdir()
        print(f'START L={box:g}, N={count}',flush=True)
        subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker','--box',str(box),'--ntracers',str(count),
                        '--output',str(dest),'--paper-source',str(args.paper_source)],check=True)
        results.append(json.loads((dest/'summary.json').read_text()))
        (root/'partial.json').write_text(json.dumps(results,indent=2)+'\n')
    sources=[Path(__file__),Path(__file__).with_name('spatial_engulfment.py'),Path(__file__).with_name('benchmark_tracers.py'),
             Path(__file__).with_name('tracer_reference.py'),Path(__file__).with_name('OPTIMIZATION_EXPERIMENT.md'),
             Path(__file__).with_name('analytic.py'),Path(__file__).with_name('run_survival.py'),ROOT/'src/axion_fopt_next/core.py']
    result=dict(status='passed',configurations=results,physics=dict(h=1.,beta=1.,speed=.5,end=2.5,dt=.05),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                limits='Geometry only; same float32 semantics as pinned reference, including existing finite-tracer/timestep bias. Three paired seeds per configuration.')
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in result['source_sha256']:
        target=root/'sources'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)


if __name__=='__main__':main()
