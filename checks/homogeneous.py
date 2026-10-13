"""Controlled radiation-era homogeneous lattice checks; run from repo root."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.special import gamma, jv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from axion_fopt_next.interchange import write_arrival_map

BOX = 2 * np.pi
FINAL = 8.


def bessel(u, amplitude):
    """m=H*=omega=1: t=tau^2/2, a=tau=1+u, pi=a^3*theta_t."""
    a = 1 + np.asarray(u)
    z = a**2 / 2
    factor = amplitude * 2**.25 * gamma(1.25) * z**(-.25)
    return factor * jv(.25, z), -a**3 * factor * jv(1.25, z)


def ode(u, amplitude, linear=False):
    initial = np.array(bessel(0., amplitude))
    def rhs(time, state):
        a = 1 + time
        theta, pi = state
        return [pi/a**2, -a**4 * (theta if linear else np.sin(theta))]
    solution = solve_ivp(rhs, (0., FINAL), initial, method="DOP853",
                         rtol=2e-12, atol=amplitude*2e-14,
                         max_step=.01, dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.sol(u)


def error(actual, expected, amplitude):
    return float(np.max(np.abs(actual - expected)) / amplitude)


def run(executable, output, n, dt, amplitude):
    name = f"N{n}_dt{dt:g}_theta{amplitude:g}"
    path = output / name
    arrival = np.zeros((n, n, n))
    map_path = output / f"{name}.afopt"
    write_arrival_map(map_path, arrival, BOX)
    theta, pi = bessel(0., amplitude)
    params = dict(N=n, lSide=BOX, dt=dt, tMax=FINAL, evolver="PV2",
                  expansion="true", fixedBackground="true", H0=1.,
                  omegaEoS=1/3, tau_origin=1., omega_star=1., f_axion=1.,
                  theta_initial=float(theta), theta_prime_initial=float(pi),
                  mass_false=0., mass_true=1., arrival_map=str(map_path),
                  result_dir=str(path), output_every=1)
    (output / f"{name}.input.json").write_text(json.dumps(params, indent=2)+"\n")
    print(f"Running {name}", flush=True)
    result = subprocess.run([str(executable), *[f"{k}={v}" for k, v in params.items()]],
                            capture_output=True, text=True, timeout=180)
    (output / f"{name}.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{name} failed; see {output / (name+'.log')}")
    history = np.genfromtxt(path / "history.csv", delimiter=",", names=True)
    initial = np.genfromtxt(path / "initial.csv", delimiter=",", names=True)
    final = np.genfromtxt(path / "final.csv", delimiter=",", names=True)
    expected_indices = np.indices((n,n,n)).reshape(3,-1).T
    for snapshot in (initial, final):
        np.testing.assert_array_equal(np.column_stack([snapshot[c] for c in ("x","y","z")]), expected_indices)
        np.testing.assert_array_equal(snapshot["arrival"], 0.)
        np.testing.assert_array_equal(snapshot["mass2"], 1.)
    np.testing.assert_allclose(initial["theta"], theta, rtol=1e-13)
    np.testing.assert_allclose(initial["pi"], pi, rtol=1e-13)
    u = history["u"]
    np.testing.assert_allclose(u, np.linspace(0., FINAL, round(FINAL/dt)+1), atol=1e-12, rtol=0)
    np.testing.assert_allclose(history["tau"], 1+u, atol=1e-12, rtol=0)
    np.testing.assert_allclose(history["a"], 1+u, atol=1e-12, rtol=0)
    analytic = np.array(bessel(u, amplitude))
    nonlinear = ode(u, amplitude)
    linear = ode(u, amplitude, linear=True)
    np.savetxt(path / "reference.csv", np.column_stack([u, (1+u)**2/2, *analytic, *nonlinear, *linear]),
               delimiter=",", header="u,mt,bessel_theta,bessel_pi,cosine_theta,cosine_pi,linear_ode_theta,linear_ode_pi", comments="")
    record = dict(name=name, n=n, dt=dt, amplitude=amplitude,
                  theta_error_cosine=error(history["theta_mean"], nonlinear[0], amplitude),
                  pi_error_cosine=error(history["pi_mean"], nonlinear[1], amplitude),
                  theta_error_bessel=error(history["theta_mean"], analytic[0], amplitude),
                  pi_error_bessel=error(history["pi_mean"], analytic[1], amplitude),
                  nonlinear_theta_correction=error(nonlinear[0], analytic[0], amplitude),
                  linear_theta_error=error(linear[0], analytic[0], amplitude),
                  linear_pi_error=error(linear[1], analytic[1], amplitude),
                  local_theta_spread=max(float(np.ptp(s["theta"]))/amplitude for s in (initial, final)),
                  local_pi_spread=max(float(np.ptp(s["pi"]))/amplitude for s in (initial, final)),
                  max_moment_variance_normalized=float(np.max(np.abs(history["theta2_mean"]-history["theta_mean"]**2))/amplitude**2))
    return record, history, analytic, nonlinear


def plots(output, runs):
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/afopt-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Scientific diagnostic plots: time-series/refinement templates with explicit
    # colors, line styles, fixed normalization, and full sampled reference CSVs.
    plt.rcParams.update({"font.size":10, "axes.spines.top":False,
                         "axes.spines.right":False, "axes.grid":True,
                         "grid.color":"#E6E8F0", "grid.linewidth":.8})
    grids = runs[:4]
    fig, axes = plt.subplots(3,1,figsize=(10,9),sharex=True)
    fig.subplots_adjust(top=.86, hspace=.3, bottom=.09)
    fig.text(.125,.97,"Homogeneous cosine evolution and grid independence",fontsize=14,weight="bold")
    fig.text(.125,.935,"Radiation era · m = H* = 1 · L = 2π/m · θᵢ = 0.01 · Δτ = 0.002",fontsize=10)
    r, h, analytic, nonlinear = grids[0]
    time = (1+h["u"])**2/2
    axes[0].plot(time, analytic[0]/r["amplitude"], color="#1F2430",lw=1.5,label="Linear Bessel")
    axes[0].plot(time, nonlinear[0]/r["amplitude"], color="#F0986E",ls="--",lw=1.5,label="Cosine ODE")
    styles = ["o", "s", "^", "x"]
    for (rec, hist, ana, ode_ref), marker in zip(grids, styles):
        axes[0].plot(time[::100], hist["theta_mean"][::100]/rec["amplitude"],
                     color="#5477C4",marker=marker,ms=4,fillstyle="none",ls="none",label=f"N = {rec['n']}")
        axes[1].plot(time, (hist["theta_mean"]-ana[0])/rec["amplitude"],
                     color="#5477C4",ls=["-","--",":","-."][styles.index(marker)],label=f"N = {rec['n']}")
        axes[2].plot(time, (hist["theta_mean"]-ode_ref[0])/rec["amplitude"],
                     color="#5477C4",ls=["-","--",":","-."][styles.index(marker)])
    axes[0].legend(ncol=3,frameon=False,fontsize=9)
    axes[0].set_ylabel(r"$\langle\theta\rangle/\theta_i$")
    axes[1].set_ylabel(r"$(\theta_{lat}-\theta_B)/\theta_i$")
    axes[2].set_ylabel(r"$(\theta_{lat}-\theta_{ODE})/\theta_i$")
    axes[2].set_xlabel(r"Cosmic time $mt$")
    for ax in axes:
        ax.set_xlim(time[0],time[-1]); ax.axhline(0,color="#6F768A",lw=.6)
    for suffix in ("png","pdf"):
        fig.savefig(output / f"homogeneous.{suffix}",dpi=180)
    plt.close(fig)
    refined = [runs[1], runs[4], runs[5]]
    dts = np.array([item[0]["dt"] for item in refined])
    errors = np.array([item[0]["theta_error_cosine"] for item in refined])
    fig, ax = plt.subplots(figsize=(7.5,5))
    fig.subplots_adjust(top=.8,bottom=.14,left=.15)
    fig.text(.15,.96,"Temporal refinement of the homogeneous mode",fontsize=13,weight="bold")
    fig.text(.15,.91,"N = 8 · θᵢ = 0.01 · max field error normalized by θᵢ · cosine ODE reference",fontsize=9)
    ax.loglog(dts,errors,"o-",color="#5477C4",label="Measured error")
    ax.loglog(dts,errors[0]*(dts/dts[0])**2,"--",color="#1F2430",label="Δτ² guide")
    ax.set_xlabel(r"Conformal timestep $\Delta\tau$"); ax.set_ylabel(r"$\max|\theta_{lat}-\theta_{ODE}|/\theta_i$")
    ax.legend(frameon=False)
    for suffix in ("png","pdf"):
        fig.savefig(output / f"temporal_refinement.{suffix}",dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable",type=Path,default=ROOT / "build/cosmolattice/axion_fopt_evolve")
    parser.add_argument("--output",type=Path,default=ROOT / "outputs/homogeneous-checks")
    options = parser.parse_args()
    executable, output = options.executable.resolve(), options.output.resolve()
    if not executable.is_file():
        parser.error(f"Missing executable: {executable}")
    output.mkdir(parents=True,exist_ok=False)
    configs = [(n,.002,.01) for n in (4,8,16,32)] + [(8,.001,.01),(8,.0005,.01),(8,.0005,.001)]
    runs = [run(executable,output,*config) for config in configs]
    gates = []
    def gate(name, measured, limit, greater=False):
        passed = bool(np.isfinite(measured) and (measured > limit if greater else measured < limit))
        gates.append(dict(name=name, measured=measured, operator=">" if greater else "<", limit=limit, passed=passed))
    for rec, *_ in runs:
        for key in ("linear_theta_error","linear_pi_error"):
            gate(f"{rec['name']}/{key}",rec[key],1e-8)
        gate(f"{rec['name']}/cosine_field",rec["theta_error_cosine"],2e-5 if rec["dt"]==.0005 else 3e-4)
        gate(f"{rec['name']}/bessel_field",rec["theta_error_bessel"],5e-4)
        for key in ("local_theta_spread","local_pi_spread"):
            gate(f"{rec['name']}/{key}",rec[key],1e-10)
        gate(f"{rec['name']}/moment_variance",rec["max_moment_variance_normalized"],1e-10)
    baseline = runs[0][1]
    for rec, history, *_ in runs[:4]:
        for field in ("theta_mean","pi_mean"):
            gate(f"N{rec['n']}/grid_{field}",error(history[field],baseline[field],.01),1e-10)
    ratio = runs[1][0]["theta_error_cosine"]/runs[5][0]["theta_error_cosine"]
    gate("temporal_coarse_to_fine_ratio",ratio,10.,greater=True)
    ratio = runs[5][0]["nonlinear_theta_correction"]/runs[6][0]["nonlinear_theta_correction"]
    gate("nonlinear_amplitude_ratio",ratio,80.,greater=True)
    files = [Path(__file__), ROOT/"src/axion_fopt_next/core.py", ROOT/"src/axion_fopt_next/interchange.py",
             *[ROOT/"adapters/cosmolattice"/name for name in ("evolve.cpp","axion_fopt.h","potential_functions.h","arrival_map.hpp","prepare_overlay.py","CMakeLists.txt")]]
    summary = dict(status="passed" if all(g["passed"] for g in gates) else "failed",
                   source="https://arxiv.org/pdf/2607.01333, Appendix B, B1-B2",
                   conventions=dict(m=1,H_star=1,L=BOX,tau_initial=1,t_initial=.5,u_final=FINAL,t_final=40.5,
                                    initial_conditions="Bessel field and derivative at finite start; theta_i is early-time regular amplitude"),
                   normalization="max absolute error / fixed theta_i; pi=a^2 dtheta/du",
                   python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,
                   executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                   source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                   git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                   git_status=subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True),
                   runs=[r[0] for r in runs],gates=gates,
                   limits="Homogeneous mode only; no bubble, non-expanding run, spatial convergence or paper reproduction claim")
    build = executable.parent/"build-info.json"
    if build.exists(): summary["build"] = json.loads(build.read_text())
    (output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    plots(output,runs)
    print(json.dumps({"status":summary["status"],"failed_gates":[g for g in gates if not g["passed"]],"output":str(output)},indent=2))
    return 0 if summary["status"]=="passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
