"""Hubble-rate survival scan; parameters/gates in HUBBLE_EXPERIMENT.md."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import time
import numpy as np
from scipy.integrate import trapezoid
from run_survival import ROOT, curve, statistics, noncausal_births
from axion_fopt_next import sample_candidates, exclude_converted_births
from analytic import Case, survival, cosmic_survival, tabulated_survival
from tracer_reference import DEFAULT_SOURCE, load_paper, simulate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/stochastic/hubble-reference')
    parser.add_argument('--paper-source',type=Path,default=DEFAULT_SOURCE)
    args=parser.parse_args()
    root=args.output.resolve(); root.mkdir(parents=True,exist_ok=False)
    shutil.copy2(args.paper_source,root/'paper_axion_sim.py.txt')
    paper,provenance=load_paper(args.paper_source)
    import numba
    numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
    cases=[Case(f'ratio1_H{h:g}',3/h,.5,2.5/h,h=h,beta=h) for h in (.5,1.,2.)]
    cases += [Case(f'beta1_H{h:g}',3/h,.5,2.5/h,h=h,beta=1.) for h in (.75,1.5)]
    gates=[]; reports=[]; scaled_arrivals={}; scaled_reference=None
    def gate(name,value,limit):
        gates.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(np.isfinite(value) and value<=limit)))
    for case in cases:
        tau=np.linspace(0,case.end,26); x=case.h*case.cosmic_time(tau)
        physical=survival(case,tau)
        gate(case.name+'/cosmic_conformal',np.max(np.abs(physical-cosmic_survival(case,tau))),1e-9)
        nodes=np.linspace(0,case.end,257); q=case.rate(nodes)
        expected=case.box**3*trapezoid(q,nodes)
        if expected>100_000: raise ValueError('Candidate budget exceeded')
        if case.speed*case.end>=case.box/2: raise ValueError('Euclidean sphere reference invalid')
        linear=tabulated_survival(nodes,q,tau,case.speed)
        fine_nodes=np.linspace(0,case.end,513)
        finer=tabulated_survival(fine_nodes,case.rate(fine_nodes),tau,case.speed)
        gate(case.name+'/rate_physical',np.max(np.abs(linear-physical)),.002)
        gate(case.name+'/rate_refinement',np.max(np.abs(linear-finer)),.002)
        if case.beta==case.h:
            if scaled_reference is None: scaled_reference=physical
            gate(case.name+'/scale_reference',np.max(np.abs(physical-scaled_reference)),1e-10)
        folder=root/case.name; folder.mkdir()
        np.savetxt(folder/'rate.csv',np.column_stack([nodes,q]),delimiter=',',header='tau,q',comments='')
        engines={}
        for method in ('poisson','tracer'):
            start=time.monotonic(); geometric=[]; internal=[]; halves=[]; counts=[]
            dest=folder/method; dest.mkdir()
            print(f'START {case.name}/{method}; expected candidates={expected:.1f}',flush=True)
            for seed in range(24):
                probes=np.random.default_rng(np.random.SeedSequence([92842,seed])).uniform(0,case.box,(2048,3))
                if method=='poisson':
                    candidates=sample_candidates(nodes,q,case.box,case.speed,seed,max_candidates=100_000)
                    catalog=exclude_converted_births(candidates)
                    gate(f'{case.name}/{seed}/causal',noncausal_births(catalog),0)
                    arrival=catalog.arrival_times(probes)
                    gate(f'{case.name}/{seed}/phantom',case.h*np.max(np.abs(arrival-candidates.arrival_times(probes))),1e-10)
                    if case.beta==case.h:
                        if seed not in scaled_arrivals: scaled_arrivals[seed]=case.h*arrival
                        gate(f'{case.name}/{seed}/scale_arrivals',np.max(np.abs(case.h*arrival-scaled_arrivals[seed])),1e-10)
                    inside=np.empty(0)
                else:
                    catalog,inside=simulate(paper,case,seed,2048,.05/case.h,tau)
                    arrival=catalog.arrival_times(probes); internal.append(inside)
                observed=curve(arrival,tau); geometric.append(observed)
                halves.append(curve(arrival[:1024],tau)-curve(arrival[1024:],tau))
                gate(f'{case.name}/{method}/{seed}/monotonic',max(0.,np.max(np.diff(observed))),0)
                gate(f'{case.name}/{method}/{seed}/bounds',max(0.,-observed.min(),observed.max()-1),0)
                counts.append(dict(seed=seed,candidates_or_tracer_births=catalog.candidate_count,accepted=len(catalog.times)))
                np.savez_compressed(dest/f'seed{seed}.npz',times=catalog.times,positions=catalog.positions,
                                    geometric_F=observed,internal_F=inside,probe_half_difference=halves[-1])
            stats=statistics(geometric); hs=statistics(halves)
            gate(case.name+'/'+method+'/probe_halves',np.max(np.abs(hs['mean'])-4*np.array(hs['stderr'])-.005),0)
            if method=='poisson':
                gate(case.name+'/poisson_agreement',np.max(np.abs(np.array(stats['mean'])-linear)-4*np.array(stats['stderr'])-.005),0)
            engines[method]=dict(geometric_F=stats,internal_F=statistics(internal) if internal else None,
                                 probe_half_difference=hs,counts=counts,wall_seconds=time.monotonic()-start,
                                 max_abs_mean_error=float(np.max(np.abs(np.array(stats['mean'])-physical))))
            print(f'DONE {case.name}/{method}: {engines[method]["wall_seconds"]:.2f}s; F_end={stats["mean"][-1]:.6f}',flush=True)
        reports.append(dict(case=asdict(case),tau=tau.tolist(),time_cosmic=case.cosmic_time(tau).tolist(),
                            Ht=x.tolist(),analytic_F=physical.tolist(),linear_reference_F=linear.tolist(),
                            expected_candidates=expected,engines=engines))
        (root/'partial.json').write_text(json.dumps(reports,indent=2)+'\n')
    files=[Path(__file__),Path(__file__).with_name('analytic.py'),Path(__file__).with_name('tracer_reference.py'),
           Path(__file__).with_name('run_survival.py'),ROOT/'src/axion_fopt_next/core.py',Path(__file__).with_name('HUBBLE_EXPERIMENT.md')]
    result=dict(status='passed' if all(g['passed'] for g in gates) else 'failed',cases=reports,gates=gates,
                seeds=list(range(24)),probes=2048,rate_nodes=257,ntracers=2048,dimensionless_dt=.05,
                paper_source=provenance,numpy=np.__version__,numba_threads=numba.get_num_threads(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(root,reports)
    print(json.dumps(dict(status=result['status'],failed=[g for g in gates if not g['passed']]),indent=2))


def plot(root,reports):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True,sharey='row',layout='constrained')
    fig.suptitle('Hubble-rate checks: geometric false-vacuum fraction\n24 realizations · bands show ±1 ensemble standard error',fontsize=14)
    for row,group in enumerate((reports[:3],[reports[3],reports[1],reports[4]])):
        for ax,r in zip(axes[row],group):
            c=r['case']; x=np.array(r['Ht']); ref=np.array(r['analytic_F'])
            ax.plot(x,ref,color='black',label='Analytic')
            for method,color in (('poisson','#4477aa'),('tracer','#cc6677')):
                stats=r['engines'][method]['geometric_F']; m=np.array(stats['mean']); se=np.array(stats['stderr'])
                ax.plot(x,m,color=color,ls='--',label=method.capitalize()+' geometry')
                ax.fill_between(x,m-se,m+se,color=color,alpha=.18)
            internal=r['engines']['tracer']['internal_F']['mean']
            ax.plot(x,internal,color='#999933',ls=':',label='Tracer internal')
            ax.set_title(f'H*={c["h"]:g}, beta={c["beta"]:g}, beta/H*={c["beta"]/c["h"]:.3g}')
            ax.set_ylim(-.03,1.03); ax.grid(alpha=.2)
            if row==1: ax.set_xlabel('H* t (cosmic time)')
        axes[row,0].set_ylabel(('Fixed beta/H*=1' if row==0 else 'Fixed beta=1')+'\nFalse-vacuum fraction F')
    axes[0,0].legend(fontsize=8)
    for ext in ('png','pdf'): fig.savefig(root/f'hubble_survival.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__': main()
