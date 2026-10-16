"""Predeclared correctness gates for the compiled dilaton+axion pocket model."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from model import (EXECUTABLE, config, run, read_state, read_slice, history, reference, potential,
                   force, dilaton_potential, lambda_for_width, tension, energies, gas_modes, gas_field, grid)

GATES = dict(force_difference=2e-7, vacuum_mass_error=2e-6, kink_tension_relative=1e-8,
             kink_residual=1e-12, gas_direct_sum=1e-10, final_ode_error=2e-4, refinement_ratio=3.5,
             snapshot_energy_relative=1e-12, energy_drift_relative=1e-5)
LAM, DV, MASS = 2., .05, 2.


def small_state(n, box):
    x = grid(n, box);X, Y, Z = np.meshgrid(x, x, x, indexing='ij');k = 2*np.pi/box
    s = np.zeros((n, n, n, 4))
    s[..., 0] = .5+.4*np.sin(k*X)+.1*np.cos(k*Y)
    s[..., 1] = .3*np.cos(k*Z)+.2*np.sin(k*(X+Y))
    s[..., 2] = .05*np.cos(k*Y)
    s[..., 3] = -.1*np.sin(k*Z)
    return s


def main():
    p = argparse.ArgumentParser();p.add_argument('--output', type=Path, required=True)
    p.add_argument('--executable', type=Path, default=EXECUTABLE);args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'contract.json').write_text(json.dumps(GATES, indent=2)+'\n')
    checks, details = {}, {}
    rng = np.random.default_rng(0);phi = rng.uniform(-.5, 1.5, 200);a = rng.normal(size=200);eps = 1e-6
    numerical = ((potential(phi+eps, a, LAM, DV, MASS**2)-potential(phi-eps, a, LAM, DV, MASS**2))/(2*eps),
                 (potential(phi, a+eps, LAM, DV, MASS**2)-potential(phi, a-eps, LAM, DV, MASS**2))/(2*eps))
    analytic = force(phi, a, LAM, DV, MASS**2)
    err = float(max(np.max(np.abs(n-f)) for n, f in zip(numerical, analytic)))
    details['force_difference'] = err;checks['force_difference'] = err < GATES['force_difference']
    masses = []
    for vacuum, expected in ((0., (LAM/2, 0.)), (1., (LAM/2+12*DV, MASS**2))):
        h = [(force(vacuum+eps, 0., LAM, DV, MASS**2)[0]-force(vacuum-eps, 0., LAM, DV, MASS**2)[0])/(2*eps),
             (force(vacuum, eps, LAM, DV, MASS**2)[1]-force(vacuum, -eps, LAM, DV, MASS**2)[1])/(2*eps)]
        masses.append(max(abs(x-y) for x, y in zip(h, expected)))
    details['vacuum_mass_error'] = max(masses);checks['vacuum_masses'] = max(masses) < GATES['vacuum_mass_error']
    checks['vacuum_energies'] = bool(dilaton_potential(1., LAM, DV) == 0. and abs(dilaton_potential(0., LAM, DV)-DV) < 1e-15)
    # Degenerate (dv=0) kink: static residual and tension used by the thin-wall references.
    width = 2.;lam = lambda_for_width(width);z = np.linspace(-40, 40, 400001)
    t = np.tanh(z/width);sech2 = 1-t*t;kink = .5*(1+t)
    residual = -sech2*t/width**2-force(kink, 0., lam, 0., 0.)[0]
    energy_density = .5*(sech2/(2*width))**2+dilaton_potential(kink, lam, 0.)
    numeric_tension = np.trapezoid(energy_density, z)
    details['kink_residual'] = float(np.max(np.abs(residual)))
    details['kink_tension_relative'] = abs(numeric_tension/tension(lam)-1)
    checks['kink_residual'] = details['kink_residual'] < GATES['kink_residual']
    checks['kink_tension'] = details['kink_tension_relative'] < GATES['kink_tension_relative']
    # Gas field equals a direct mode sum at cell centres.
    modes = gas_modes(16., .5, 1.3, 3);n = 24;a_grid, p_grid = gas_field(modes, n);x = grid(n, 16.)
    k = 2*np.pi*modes['vectors']/16.;picks = rng.integers(0, n, (20, 3));direct = []
    for i, j, l in picks:
        phase = np.exp(1j*(k@np.array([x[i], x[j], x[l]])))
        direct.append(max(abs((modes['field']*phase).sum().real-a_grid[i, j, l]),
                          abs((modes['momentum']*phase).sum().real-p_grid[i, j, l])))
    details['gas_direct_sum'] = float(max(direct));checks['gas_direct_sum'] = max(direct) < GATES['gas_direct_sum']
    # Compiled PV2 versus independent semi-discrete reference.
    n, box = 6, 6.;s = small_state(n, box);errors = []
    for dt in (.005, .0025):
        cfg = config(n, box, dt, 1., lam=LAM, dv=DV, mass_tv=MASS, every=1, slices=int(round(1./dt)))
        directory = args.output/f'small_{dt}';cfg = run(s, cfg, directory, args.executable)
        h = history(directory);ref = reference(s, cfg, [0., 1.])[1]
        meta, final = read_state(directory/'native/final.bin');_, initial = read_state(directory/'native/initial.bin')
        error = float(np.max(np.abs(final-ref)));errors.append(error)
        checks[f'small_{dt}_transfer'] = bool(np.array_equal(initial, s))
        checks[f'small_{dt}_time'] = abs(meta['u']-1.) < 1e-14
        independent = energies(final, box, LAM, DV, MASS**2)
        scale = max(abs(independent['total']), 1e-12)
        energy_err = max(abs(independent[key]-h[key][-1]) for key in independent if key in h.dtype.names)/scale
        checks[f'small_{dt}_energy'] = energy_err < GATES['snapshot_energy_relative']
        drift = float(np.max(np.abs(h['total']-h['total'][0]))/abs(h['total'][0]))
        checks[f'small_{dt}_conservation'] = drift < GATES['energy_drift_relative']
        _, first = read_slice(directory/'native/slice_0000000.bin')
        _, last = read_slice(directory/f'native/slice_{int(round(1./dt)):07d}.bin')
        checks[f'small_{dt}_slices'] = bool(np.array_equal(first, initial[:, :, n//2]) and np.array_equal(last, final[:, :, n//2]))
        details[f'small_{dt}'] = dict(final_ode_error=error, energy_error=energy_err, energy_drift=drift)
    checks['small_ode'] = max(errors) < GATES['final_ode_error']
    details['refinement_ratio'] = errors[0]/errors[1]
    checks['small_refinement'] = errors[0]/errors[1] > GATES['refinement_ratio']
    # Bad input must fail without creating native output.
    bad = args.output/'bad.in';base = (args.output/'small_0.005/input.in').read_text()
    for label, replacement in [('nonfinite', 'lambda = nan'), ('fractional_grid', 'N = 6.5'),
                               ('expansion', 'expansion = true'), ('negative_mass', 'mass_tv = -1')]:
        key = replacement.split(' = ')[0]
        text = '\n'.join(replacement if line.startswith(key+' = ') else line for line in base.splitlines())
        text = text.replace(str((args.output/'small_0.005/native').resolve()), str((args.output/label).resolve()))
        bad.write_text(text+'\n')
        result = subprocess.run([str(args.executable), 'input='+str(bad.resolve())], capture_output=True, text=True)
        checks['reject_'+label] = result.returncode != 0 and not (args.output/label).exists()
    checks = {k: bool(v) for k, v in checks.items()}
    report = dict(gates=GATES, checks=checks, details=details, passed=all(checks.values()),
                  executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest())
    (args.output/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2));raise SystemExit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
