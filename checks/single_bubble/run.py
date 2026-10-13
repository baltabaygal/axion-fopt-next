"""Reproduce the smallest Table-I boxes of arXiv:2402.09501, Sec. IV F."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import scipy
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from axion_fopt_next import Catalogue
from axion_fopt_next.interchange import write_arrival_map
from reference import mass, observables, read_state

BOX = 2*np.pi


def case(executable,root,n,v,refinement=1):
    name=f"N{n}_v{v:g}_ref{refinement}"
    end=BOX/v
    # Fixed dt across N for each v; finest N=128 sets the stability bound.
    steps=int(np.ceil(end/(.2*BOX/128)))*refinement
    dt=end/steps
    arrival=Catalogue(BOX,v,[0.],[[BOX/2]*3]).lattice_arrival_times(n)
    map_path=root/f"{name}.afopt"
    write_arrival_map(map_path,arrival,BOX)
    out=root/name
    params=dict(N=n,lSide=BOX,dt=dt,tMax=end,evolver="PV2",expansion="false",
                theta_initial=1.,theta_prime_initial=0.,mass_false=0.,mass_true=1.,
                potential="quadratic",wall_width=.15,profile_wall_speed=v,
                omega_star=1.,f_axion=1.,tau_origin=0.,arrival_map=str(map_path),
                result_dir=str(out),snapshot_format="binary",energy_diagnostics="true",
                output_every=max(1,round(.02/dt)))
    if n==128 and v==.5 and refinement==1:
        params["snapshot_at"]=2.
    (root/f"{name}.input.json").write_text(json.dumps(params,indent=2)+"\n")
    print(f"START {name}: {steps} steps, dt={dt:.6g}",flush=True)
    begin=time.monotonic()
    with (root/f"{name}.log").open("w") as log:
        result=subprocess.run([str(executable),*[f"{k}={val}" for k,val in params.items()]],stdout=log,stderr=subprocess.STDOUT,timeout=7200)
    if result.returncode:
        raise RuntimeError(f"{name} failed: see log")
    meta,initial=read_state(out/"initial.bin")
    final_meta,final=read_state(out/"final.bin")
    np.testing.assert_array_equal(initial[...,0],1.)
    np.testing.assert_array_equal(initial[...,1],0.)
    np.testing.assert_array_equal(initial[...,2],arrival)
    for t,state in ((0.,initial),(end,final)):
        np.testing.assert_allclose(state[...,3],mass(t,arrival,v),atol=2e-14,rtol=0.)
    np.testing.assert_allclose(final_meta["time"],end,rtol=0,atol=1e-12)
    if meta["scale"]!=1 or final_meta["scale"]!=1:
        raise RuntimeError("Unexpected expansion")
    history=np.genfromtxt(out/"history.csv",delimiter=",",names=True)
    np.testing.assert_allclose(history["a"],1.,atol=1e-14)
    stats,spectrum=observables(final[...,0],final[...,1],BOX)
    np.savetxt(out/"spectrum.csv",spectrum,delimiter=",",header="k_over_m0,dn_dlogk_over_m0_phi0_squared",comments="")
    injection=trapezoid(history["source_mean"],history["u"])
    balance=history["total_mean"][-1]-history["total_mean"][0]-injection
    fit=(1.3*v**.4 if v<=.4 else 1.8*v**.8)/BOX
    record=dict(name=name,n=n,v=v,refinement=refinement,dt=dt,steps=steps,end=end,
                wall_width=.15,wall_cells=.15/(BOX/n),elapsed_seconds=time.monotonic()-begin,
                **stats,final_energy=float(history["total_mean"][-1]),
                energy_work_relative_residual=float(abs(balance)/history["total_mean"][-1]),
                minimum_final_mass_squared=float(final[...,3].min()),
                paper_piecewise_fit=fit,relative_difference_from_fit=float((stats["number_continuum"]-fit)/fit))
    (out/"observables.json").write_text(json.dumps(record,indent=2)+"\n")
    print(f"DONE {name}: n/(m0 phi0²)={stats['number_continuum']:.8f}; {record['elapsed_seconds']:.1f}s",flush=True)
    return record


def summarize(root,records,executable):
    gates=[]
    def gate(name,value,limit):
        gates.append(dict(name=name,value=float(value),limit=limit,passed=bool(np.isfinite(value) and value<limit)))
    for r in records:
        gate(r["name"]+"/parseval",r["parseval_relative_error"],1e-12)
        gate(r["name"]+"/external_work_balance",r["energy_work_relative_residual"],.005)
        gate(r["name"]+"/final_mass",abs(1-r["minimum_final_mass_squared"]),1e-4)
    indexed={(r["n"],r["v"],r["refinement"]):r for r in records}
    for n in (64,128):
        if (n,.5,1) in indexed and (n,.5,2) in indexed:
            coarse=indexed[n,.5,1]; fine=indexed[n,.5,2]
            gate(f"N{n}/time_number_change",abs(coarse["number_continuum"]/fine["number_continuum"]-1),.002)
    if (64,.5,1) in indexed and (128,.5,1) in indexed:
        gate("v0.5/space_number_change",abs(indexed[64,.5,1]["number_continuum"]/indexed[128,.5,1]["number_continuum"]-1),.02)
    sources=[Path(__file__),Path(__file__).with_name("reference.py"),ROOT/"src/axion_fopt_next/core.py",
             *[ROOT/"adapters/cosmolattice"/p for p in ("axion_fopt.h","potential_functions.h","evolve.cpp","prepare_overlay.py","CMakeLists.txt")]]
    result=dict(status="passed" if all(g["passed"] for g in gates) else "failed",runs=records,gates=gates,
                paper="https://arxiv.org/pdf/2402.09501 (v2), Sec. IV F, Eqs.53-57, Table I",
                conventions="m0=phi0=omega_star=f_axion=1; a=1; t=u=tau; box L=2pi; center=L/2; initial phi=1, phi_dot=0; delta_w=.15",
                fit_caveat="Eqs.56-57 approximate fit measured at L=16pi; used as a comparison, not an exact acceptance target for L=2pi",
                limits="Only smallest Table-I box; PV2 and nearest-neighbor stencil; cell-centered sampling; no digitized figure or author raw data; spatial/time refinements at v=.5 only",
                executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                numpy=np.__version__,scipy=scipy.__version__,python=sys.version.split()[0],
                git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                git_status=subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True),
                build=json.loads((executable.parent/"build-info.json").read_text()))
    (root/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


def plots(root,records):
    os.environ.setdefault("MPLCONFIGDIR","/tmp/afopt-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.alpha":.2})
    paper=sorted([r for r in records if r["n"]==128 and r["refinement"]==1],key=lambda r:r["v"])
    if not paper: return
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    fig.subplots_adjust(top=.8,bottom=.16,wspace=.3)
    fig.text(.125,.96,"Centered periodic bubble: non-expanding harmonic field",fontsize=14,weight="bold")
    fig.text(.125,.9,"L = 2π/m₀ · N = 128 · δw = 0.15/m₀ · φ(t=0) = φ₀ · φ̇(t=0) = 0",fontsize=10)
    for r,color,style in zip(paper,("#1F2430","#5477C4","#CC6F47"),("-","--","-.")):
        h=np.genfromtxt(root/r["name"]/"history.csv",delimiter=",",names=True)
        axes[0].plot(h["u"]/r["end"],h["total_mean"],color=color,ls=style,label=f"v = {r['v']}")
        axes[0].plot(h["u"]/r["end"],h["inside_mean"],color=color,ls=":",alpha=.75)
        axes[1].plot(r["v"],r["number_continuum"],"o",color=color,ms=7)
    vel=np.linspace(.1,.9,150)
    fit=np.where(vel<=.4,1.3*vel**.4,1.8*vel**.8)/BOX
    axes[1].plot(vel,fit,color="#6F768A",ls="--",label="Paper fit (L = 16π/m₀)")
    axes[0].set_xlabel(r"$t/t_{end}$"); axes[0].set_ylabel(r"$\bar\rho/(m_0^2\phi_0^2)$")
    axes[0].legend(frameon=False); axes[0].text(.04,.94,"Dotted: energy inside r < vt / box volume",transform=axes[0].transAxes,fontsize=8)
    axes[1].set_xlabel("Wall speed v"); axes[1].set_ylabel(r"$n/(m_0\phi_0^2)$"); axes[1].legend(frameon=False,fontsize=9)
    for fmt in ("png","pdf"): fig.savefig(root/f"bubble_observables.{fmt}",dpi=180)
    plt.close(fig)
    target=next((r for r in paper if r["v"]==.5),None)
    if target and (root/target["name"]/"checkpoint.bin").exists():
        meta,early=read_state(root/target["name"]/"checkpoint.bin")
        _,final=read_state(root/target["name"]/"final.bin")
        z=64
        fig,axes=plt.subplots(1,2,figsize=(10,4.8))
        fig.subplots_adjust(top=.78,bottom=.13,wspace=.35)
        fig.text(.125,.95,"Axion depletion and subsequent periodic evolution",fontsize=13,weight="bold")
        fig.text(.125,.89,"v = 0.5 · N = 128 · slice nearest z = L/2 · field normalized by φ₀",fontsize=10)
        bound=max(np.abs(early[...,0]).max(),np.abs(final[...,0]).max())
        for ax,state,title in zip(axes,(early,final),(f"mt = {meta['time']:.3f} (before overlap)",f"mt = {target['end']:.3f} (final)")):
            ax.grid(False)
            im=ax.imshow(state[:,:,z,0].T,origin="lower",extent=(0,BOX,0,BOX),cmap="RdBu_r",vmin=-bound,vmax=bound)
            ax.set_title(title,fontsize=10); ax.set_xlabel("m₀ x"); ax.set_ylabel("m₀ y")
            fig.colorbar(im,ax=ax,shrink=.85,label="φ/φ₀")
        dz=(z+.5)*BOX/128-BOX/2
        plane_radius=np.sqrt(max(0.,(.5*meta["time"])**2-dz**2))
        circle=plt.Circle((BOX/2,BOX/2),plane_radius,fill=False,color="#1F2430",ls="--",lw=1)
        axes[0].add_patch(circle)
        for fmt in ("png","pdf"): fig.savefig(root/f"bubble_slices.{fmt}",dpi=180)
        plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--executable",type=Path,default=ROOT/"build/cosmolattice/axion_fopt_evolve")
    p.add_argument("--output",type=Path,default=ROOT/"outputs/lee-single-bubble")
    p.add_argument("--suite",choices=("pilot","paper"),default="paper")
    args=p.parse_args()
    out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    exe=args.executable.resolve()
    configs=[(32,.5,1),(64,.5,1),(64,.5,2)]
    if args.suite=="paper": configs += [(128,.5,1),(128,.5,2),(128,.9,1),(128,.1,1)]
    records=[]
    for config in configs:
        records.append(case(exe,out,*config))
        summarize(out,records,exe) # preserve measured results between long runs
    summary=summarize(out,records,exe)
    plots(out,records)
    print(json.dumps({"status":summary["status"],"failed_gates":[g for g in summary["gates"] if not g["passed"]],"output":str(out)},indent=2))
    return 0 if summary["status"]=="passed" else 1


if __name__=="__main__": raise SystemExit(main())
