"""Geometry-only ensemble comparison with analytic F(t); see EXPERIMENT.md."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
import numpy as np
from scipy.integrate import trapezoid

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from axion_fopt_next import sample_candidates,exclude_converted_births
from analytic import CASES,survival,cosmic_survival,tabulated_survival
from tracer_reference import DEFAULT_SOURCE,load_paper,simulate


def curve(arrivals,times):
    return np.mean(arrivals[None,:]>times[:,None],axis=1)


def statistics(values):
    values=np.asarray(values)
    return dict(mean=values.mean(axis=0).tolist(),spread=values.std(axis=0,ddof=1).tolist(),
                stderr=(values.std(axis=0,ddof=1)/np.sqrt(len(values))).tolist())


def noncausal_births(catalogue):
    """Audit each birth against earlier periodic fronts using arrival times."""
    violations=0
    for i in range(1,len(catalogue.times)):
        delta=np.abs(catalogue.positions[:i]-catalogue.positions[i])
        delta=np.minimum(delta,catalogue.box_size-delta)
        earlier_arrival=catalogue.times[:i]+np.linalg.norm(delta,axis=1)/catalogue.wall_speed
        violations+=int(np.any(earlier_arrival<=catalogue.times[i]))
    return violations


def plot(root,reports):
    import os
    os.environ.setdefault("MPLCONFIGDIR","/tmp/afopt-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex="col")
    fig.subplots_adjust(top=.86,hspace=.28,wspace=.24,bottom=.1)
    fig.text(.125,.975,"Stochastic nucleation: geometric survival versus analytic F(t)",fontsize=14,weight="bold")
    fig.text(.125,.93,"24 realizations · independent probes · shaded band: Poisson mean ± ensemble standard error",fontsize=10)
    colors={"poisson_257":"#5477C4","tracer_512_coarse":"#CC6F47","tracer_2048_coarse":"#7A828F","tracer_2048_fine":"#1F2430"}
    styles={"poisson_257":"-","tracer_512_coarse":":","tracer_2048_coarse":"--","tracer_2048_fine":"-."}
    labels={"poisson_257":"Poisson, 257 rate nodes","tracer_512_coarse":"Tracer 512, coarse dt","tracer_2048_coarse":"Tracer 2048, coarse dt","tracer_2048_fine":"Tracer 2048, fine dt"}
    for col,report in enumerate(reports):
        case=report["case"]; x=np.array(report["time_cosmic"]); analytic=np.array(report["analytic_F"])
        axes[0,col].plot(x,analytic,color="#1F2430",lw=2,label="Analytic F")
        for label,color in colors.items():
            stats=report["engines"][label]["geometric_F"]
            mean=np.array(stats["mean"]); se=np.array(stats["stderr"])
            axes[0,col].plot(x,mean,color=color,ls=styles[label],label=labels[label])
            axes[1,col].plot(x,mean-analytic,color=color,ls=styles[label])
            if label=="poisson_257":
                axes[0,col].fill_between(x,mean-se,mean+se,color=color,alpha=.2)
                axes[1,col].fill_between(x,mean-analytic-se,mean-analytic+se,color=color,alpha=.2)
        axes[0,col].set_title(f"{case['name']} · L={case['box']:g}, v={case['speed']:g}",fontsize=11)
        axes[0,col].set_ylim(-.02,1.02); axes[0,col].set_ylabel("Geometric false-vacuum fraction F")
        axes[1,col].axhline(0,color="#6F768A",lw=.7)
        axes[1,col].set_ylabel("Ensemble mean − analytic F")
        axes[1,col].set_xlabel("Cosmic time t" if case["constant_q"] is not None else "Cosmic time H* t (H*=1)")
        for ax in axes[:,col]: ax.grid(alpha=.15); ax.set_xlim(x[0],x[-1])
    axes[0,0].legend(frameon=False,fontsize=8,loc="lower left")
    for suffix in ("png","pdf"): fig.savefig(root/f"survival.{suffix}",dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    fig.subplots_adjust(top=.8,bottom=.15,wspace=.25)
    fig.text(.125,.96,"Tracer depletion: internal estimator and independent geometry",fontsize=13,weight="bold")
    fig.text(.125,.9,"Same tracer-generated catalogues · solid: geometric F · dashed: active / initial tracers",fontsize=10)
    for ax,report in zip(axes,reports):
        x=report["time_cosmic"]
        ax.plot(x,report["analytic_F"],color="#1F2430",ls=":",label="Analytic")
        for label,color in (("tracer_512_coarse","#CC6F47"),("tracer_2048_fine","#5477C4")):
            engine=report["engines"][label]
            ax.plot(x,engine["geometric_F"]["mean"],color=color,label=f"{engine['ntracers']} tracers: geometry")
            ax.plot(x,engine["internal_F"]["mean"],color=color,ls="--",label=f"{engine['ntracers']} tracers: internal")
        ax.set_title(report["case"]["name"]); ax.set_ylim(-.02,1.02)
        ax.set_xlabel("Cosmic time t"); ax.set_ylabel("False-vacuum fraction F"); ax.grid(alpha=.15)
    axes[0].legend(frameon=False,fontsize=8)
    for suffix in ("png","pdf"): fig.savefig(root/f"tracer_estimators.{suffix}",dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"outputs/stochastic/survival-first")
    parser.add_argument("--paper-source",type=Path,default=DEFAULT_SOURCE)
    args=parser.parse_args()
    root=args.output.resolve(); root.mkdir(parents=True,exist_ok=False)
    shutil.copy2(args.paper_source,root/"paper_axion_sim.py.txt")
    paper,provenance=load_paper(args.paper_source)
    # Use the published module's thread count, without importing its side effects.
    import numba
    numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
    gates=[]; reports=[]; max_phantom_error=0.
    def gate(name,value,limit):
        gates.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(np.isfinite(value) and value<=limit)))
    for case_index,case in enumerate(CASES):
        case_root=root/case.name; case_root.mkdir()
        times=np.arange(round(case.end/.1)+1)*.1
        analytic=survival(case,times)
        cosmic=cosmic_survival(case,times)
        gate(case.name+"/analytic_time_mapping",np.max(np.abs(analytic-cosmic)),1e-9)
        if case.speed*case.end>=case.box/2: raise ValueError("Euclidean analytic volume invalid")
        engines={}
        coarse=.02 if case.constant_q is not None else .05
        configs=[("poisson_65",65,None,None),("poisson_257",257,None,None),
                 ("tracer_512_coarse",None,512,coarse),("tracer_2048_coarse",None,2048,coarse),
                 ("tracer_2048_fine",None,2048,coarse/2)]
        for label,nodes,ntracers,dt in configs:
            begin=time.monotonic(); geometric=[]; internal=[]; halves=[]; counts=[]
            catalogue_dir=case_root/label; catalogue_dir.mkdir()
            if nodes:
                tau_grid=np.linspace(0,case.end,nodes); q=case.rate(tau_grid)
                np.savetxt(catalogue_dir/"rate.csv",np.column_stack([tau_grid,q]),delimiter=",",header="tau,q",comments="")
                interpolant=tabulated_survival(tau_grid,q,times,case.speed)
                if nodes==257: gate(case.name+"/fine_rate_interpolation",np.max(np.abs(interpolant-analytic)),.002)
                if case.box**3*trapezoid(q,tau_grid)>100_000: raise ValueError("Pilot candidate budget exceeded")
            print(f"START {case.name}/{label}",flush=True)
            for seed in range(24):
                probes=np.random.default_rng(np.random.SeedSequence([92841,case_index,seed])).uniform(0,case.box,(2048,3))
                if nodes:
                    candidates=sample_candidates(tau_grid,q,case.box,case.speed,seed,max_candidates=100_000)
                    catalogue=exclude_converted_births(candidates)
                    gate(f"{case.name}/{label}/seed{seed}/causal_births",noncausal_births(catalogue),0.)
                    arrivals=catalogue.arrival_times(probes)
                    phantom=candidates.arrival_times(probes)
                    max_phantom_error=max(max_phantom_error,float(np.max(np.abs(arrivals-phantom))))
                else:
                    catalogue,tracer_curve=simulate(paper,case,seed,ntracers,dt,times)
                    arrivals=catalogue.arrival_times(probes)
                    internal.append(tracer_curve)
                full=curve(arrivals,times); half_a=curve(arrivals[:1024],times); half_b=curve(arrivals[1024:],times)
                geometric.append(full); halves.append(half_a-half_b)
                gate(f"{case.name}/{label}/seed{seed}/bounds",max(0.,-float(full.min()),float(full.max())-1.),0.)
                gate(f"{case.name}/{label}/seed{seed}/monotonic",max(0.,np.max(np.diff(full))),0.)
                if internal: gate(f"{case.name}/{label}/seed{seed}/internal_monotonic",max(0.,np.max(np.diff(internal[-1]))),0.)
                counts.append(dict(seed=seed,candidate_or_tracer_birth_count=catalogue.candidate_count,accepted_or_tracer_bubbles=len(catalogue.times)))
                np.savez_compressed(catalogue_dir/f"seed{seed}.npz",times=catalogue.times,positions=catalogue.positions,
                                    geometric_F=full,probe_half_difference=half_a-half_b,
                                    internal_F=tracer_curve if not nodes else np.empty(0))
            stats=statistics(geometric); hstats=statistics(halves)
            error=np.abs(np.array(stats["mean"])-(interpolant if nodes else analytic))
            allowance=4*np.array(stats["stderr"])+.005
            if nodes: gate(case.name+"/"+label+"/analytic_agreement",np.max(error-allowance),0.)
            gate(case.name+"/"+label+"/probe_halves",np.max(np.abs(hstats["mean"])-4*np.array(hstats["stderr"])-.005),0.)
            engines[label]=dict(ntracers=ntracers,outer_dt=dt,rate_nodes=nodes,geometric_F=stats,
                                internal_F=statistics(internal) if internal else None,probe_half_difference=hstats,
                                reference_F=interpolant.tolist() if nodes else analytic.tolist(),
                                max_abs_mean_error_physical=float(np.max(np.abs(np.array(stats["mean"])-analytic))),
                                counts=counts,wall_seconds=time.monotonic()-begin)
            print(f"DONE {case.name}/{label}: F_end={stats['mean'][-1]:.6f}, max mean error={engines[label]['max_abs_mean_error_physical']:.4f}, {engines[label]['wall_seconds']:.1f}s",flush=True)
        report=dict(case=asdict(case),tau=times.tolist(),time_cosmic=case.cosmic_time(times).tolist(),analytic_F=analytic.tolist(),engines=engines)
        reports.append(report)
        (root/"partial.json").write_text(json.dumps(reports,indent=2)+"\n")
    gate("identical_events/phantom_equivalence",max_phantom_error,1e-12)
    sources=[Path(__file__),Path(__file__).with_name("analytic.py"),Path(__file__).with_name("tracer_reference.py"),
             ROOT/"src/axion_fopt_next/core.py",Path(__file__).with_name("EXPERIMENT.md")]
    result=dict(status="passed" if all(g["passed"] for g in gates) else "failed",seeds=list(range(24)),probes=2048,
                paper_source=provenance,reference="Physical q integrated independently in conformal and cosmic time; initial cosmic time zero",
                cases=reports,gates=gates,numpy=np.__version__,numba_threads=numba.get_num_threads(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                limits="Small geometry-only ensemble; beta/H=1 radiation diagnostic; fixed imposed outer timestep; no coupled field, paper parameter scan or tracer convergence claim")
    (root/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    plot(root,reports)
    print(json.dumps({"status":result["status"],"failed_gates":[g for g in gates if not g["passed"]],"output":str(root)},indent=2))
    return 0 if result["status"]=="passed" else 1


if __name__=="__main__": raise SystemExit(main())
