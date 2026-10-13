"""Run compiled CosmoLattice against independent ODE/semi-discrete references.

Requires scipy; usage:
PYTHONPATH=src python examples/validate_cosmolattice.py --executable /path/to/axion_fopt_evolve
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import scipy
import numpy as np
from scipy.integrate import solve_ivp
from axion_fopt_next import Catalogue
from axion_fopt_next.interchange import write_arrival_map


def reference(arrival, *, dt_output, final_u, theta=.7, derivative=.13,
              origin=2., omega=2., mass_false=0., mass_true=1.3,
              expanding=False, box=4., spatial=False):
    """Integrate the canonical equations with DOP853, split at every mass jump."""
    shape = arrival.shape
    size = arrival.size if spatial else 1
    times = (arrival.ravel() if spatial else arrival.ravel()[:1])
    switches = omega * (times - origin)
    breaks = np.unique(np.r_[0., switches[(switches > 0.) & (switches < final_u)], final_u])
    y = np.r_[np.full(size, theta), np.full(size, derivative / omega)]
    samples = np.linspace(0., final_u, round(final_u / dt_output) + 1)
    result = np.empty((len(samples), 2 * size))
    for left, right in zip(breaks[:-1], breaks[1:]):
        mass2 = np.where(switches <= (left + right) / 2, mass_true**2, mass_false**2) / omega**2
        def rhs(u, state):
            a = 1. + u / (origin * omega) if expanding else 1.
            th, pi = state[:size], state[size:]
            force = -mass2 * np.sin(th)
            if spatial:
                field = th.reshape(shape)
                lap = sum(np.roll(field, 1, axis) + np.roll(field, -1, axis) - 2 * field
                          for axis in range(3)) / (omega * box / shape[0])**2
                force += lap.ravel() / a**2
            return np.r_[pi / a**2, a**4 * force]
        solution = solve_ivp(rhs, (left, right), y, method="DOP853", rtol=2e-11,
                             atol=2e-13, dense_output=True, max_step=.05)
        if not solution.success:
            raise RuntimeError(solution.message)
        select = (samples >= left) & (samples <= right)
        if np.any(select):
            result[select] = solution.sol(samples[select]).T
        y = solution.y[:, -1]
    return samples, result


def run(executable, root, name, arrival, dt, *, final=2., expanding=False,
        theta=.7, derivative=.13, mf=0., mt=1.3, omega=2., origin=2., box=4.):
    map_path = root / (name + ".afopt")
    write_arrival_map(map_path, arrival, box)
    out = root / name
    args = {
        "N": arrival.shape[0], "lSide": box * omega, "dt": dt, "tMax": final,
        "expansion": str(expanding).lower(), "fixedBackground": "true",
        "H0": 1 / origin, "omegaEoS": 1/3, "evolver": "PV2", "a0": 1,
        "theta_initial": theta, "theta_prime_initial": derivative,
        "mass_false": mf, "mass_true": mt, "omega_star": omega,
        "tau_origin": origin, "arrival_map": map_path, "result_dir": out,
    }
    result = subprocess.run([str(executable), *[f"{k}={v}" for k, v in args.items()]],
                            capture_output=True, text=True, timeout=120)
    (root / (name + ".log")).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{name} failed: {result.stderr}")
    history = np.genfromtxt(out / "history.csv", delimiter=",", names=True)
    initial = np.genfromtxt(out / "initial.csv", delimiter=",", names=True)
    end = np.genfromtxt(out / "final.csv", delimiter=",", names=True)
    np.testing.assert_array_equal(initial["arrival"], arrival.ravel())
    expected_mass = np.where(origin + final/omega >= arrival.ravel(), mt**2, mf**2)/omega**2
    np.testing.assert_allclose(end["mass2"], expected_mass, atol=1e-15)
    np.testing.assert_allclose(history["potential_mean"][-1],
                               np.mean(end["mass2"]*(1-np.cos(end["theta"]))), atol=2e-14)
    expected_a = 1 + history["u"]/(origin*omega) if expanding else np.ones(len(history))
    np.testing.assert_allclose(history["a"], expected_a, atol=2e-14)
    return history, end


def validate(executable, root):
    root.mkdir(parents=True, exist_ok=False)
    summary = {"python":sys.version.split()[0], "numpy":np.__version__, "scipy":scipy.__version__,
               "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
               "reference": "SciPy DOP853; canonical equations split at each lattice mass jump",
               "checks": {}, "limits": "CPU, single process, PV2, fixed grid; no published benchmark or spatial continuum convergence"}
    build_info = executable.parent / "build-info.json"
    if build_info.exists():
        summary["build"] = json.loads(build_info.read_text())
    empty = np.full((4,4,4), np.inf)
    # Analytic zero-mass drift, including nonzero initial momentum and RD damping.
    for expansion in (False, True):
        label = "zero_mass_rd" if expansion else "zero_mass_flat"
        history, _ = run(executable, root, label, empty, .02, expanding=expansion)
        u = history["u"]
        expected = .7 + .13 * (2 - 4/(2 + u/2)) if expansion else .7 + .13*u/2
        error = float(np.max(np.abs(history["theta_mean"] - expected)))
        assert error < 3e-6, (label, error)
        np.testing.assert_allclose(history["pi_mean"], .13/2, atol=2e-14)
        summary["checks"][label] = {"max_theta_error": error}

    # A transition before the FIRST midpoint must affect the FIRST kick.
    switch = np.full((4,4,4), 2. + .015/2)
    history, end = run(executable, root, "midpoint_switch", switch, .04,
                       final=.04, derivative=0.)
    expected_pi = -.04 * (1.3/2)**2 * np.sin(.7)
    np.testing.assert_allclose(end["pi"], expected_pi, rtol=1e-13)
    summary["checks"]["midpoint_switch"] = {"pi": float(end["pi"][0]), "expected_pi": float(expected_pi)}

    # Smooth segments; delayed homogeneous switch aligned with all step grids.
    for expansion in (False, True):
        for switch_u in (0., .4):
            arrival = np.full((4,4,4), 2. + switch_u/2)
            errors = []
            label = f"homogeneous_{'rd' if expansion else 'flat'}_switch_{switch_u}"
            for dt in (.04, .02, .01):
                history, end = run(executable, root, label + f"_dt_{dt}", arrival, dt, expanding=expansion)
                _, ref = reference(arrival, dt_output=dt, final_u=2., expanding=expansion)
                error = float(np.max(np.abs(history["theta_mean"] - ref[:,0])))
                errors.append(error)
                assert np.ptp(end["theta"]) < 1e-13, "Homogeneous solution developed spatial variation"
            assert errors[-1] < 3e-5 and errors[0]/errors[-1] > 10, (label, errors)
            summary["checks"][label] = {"dt": [.04,.02,.01], "max_theta_errors": errors,
                                        "coarse_to_fine_ratio": errors[0]/errors[-1]}

    # Non-axis-symmetric periodic bubble; compare both field and canonical momentum.
    arrival = Catalogue(4., .8, [2.1], [[.13,1.07,2.31]]).lattice_arrival_times(8)
    _, ref = reference(arrival, dt_output=.01, final_u=3., derivative=0.,
                       spatial=True, expanding=True)
    errors = []
    for dt in (.02,.01,.005):
        _, end = run(executable, root, f"one_bubble_dt_{dt}", arrival, dt,
                     final=3., derivative=0., expanding=True)
        errors.append({"dt": dt,
                       "theta_max_error": float(np.max(np.abs(end["theta"]-ref[-1,:arrival.size]))),
                       "pi_max_error": float(np.max(np.abs(end["pi"]-ref[-1,arrival.size:]))),
                       "theta_std": float(np.std(end["theta"]))})
    # Sharp wall: generic unaligned jumps reduce global time convergence order.
    assert errors[-1]["theta_max_error"] < .003 and errors[-1]["pi_max_error"] < .006, errors
    assert errors[-1]["theta_max_error"] < errors[0]["theta_max_error"], errors
    assert errors[-1]["theta_std"] > 1e-3, "One-bubble field stayed homogeneous"
    summary["checks"]["one_bubble_rd"] = {"n":8, "box_size":4., "errors":errors,
        "note":"Compared at fixed spatial resolution; sharp-wall crossings generally reduce temporal order"}
    (root / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("outputs/cosmolattice-validation"))
    options = parser.parse_args()
    validate(options.executable.resolve(), options.output.resolve())
