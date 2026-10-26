"""Field-only pressure-balanced energy pilot; contract: FIELD_GOAL.md.

The energy estimate does not make the initial random stress self-consistent.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import (EXECUTABLE, config, eigen_state, energies, history,
                   radial_modes, read_slice, read_state, run, tension)
from run_prototype import wall_cells


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    n, box, r0, width, lam, dv, mass = 128, 64., 20., 2., 2., .01, 2.
    target = 3*(dv+2*tension(lam)/r0)*4*np.pi/3*r0**3
    modes = radial_modes(r0, width, mass, .15, .6, 31., h=.01)
    s = eigen_state(n, box, r0, width, modes, 1., 0, 30.)
    e = energies(s, box, lam, dv, mass**2)
    amplitude = float(np.sqrt(target/sum(e[k] for k in ('axion_kinetic', 'axion_gradient', 'coupling'))))
    s[..., 1] *= amplitude;s[..., 3] *= amplitude
    cfg = config(n, box, .0125, 80., lam=lam, dv=dv, mass_tv=mass, every=40, slices=320)
    contract = dict(config=cfg, target_axion_energy=target, amplitude=amplitude, seed=0,
                    modes=[dict(l=m['l'], n=m['n'], omega=m['omega']) for m in modes],
                    mode_assumption='Gaussian classical high-occupation equipartition; no thermalisation',
                    sources={name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                             for name in ('FIELD_GOAL.md', 'run_balanced_pilot.py', 'model.py',
                                          'pocket_model.h', 'evolve.cpp', 'run_prototype.py')},
                    executable_sha256=hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest())
    (args.output/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    error = None
    try:
        run(s, cfg, args.output/'balanced')
    except RuntimeError as exc:
        error = str(exc)
    h = history(args.output/'balanced')
    _, initial = read_state(args.output/'balanced/native/initial.bin')
    widths = []
    for path in sorted((args.output/'balanced/native').glob('slice_*.bin')):
        meta, plane = read_slice(path)
        i = np.argmin(np.abs(h['u']-meta['u']))
        if h['fv_radius'][i] >= 8:
            widths.append(float(wall_cells(plane, box/n)))
    drift = float(np.max(np.abs(h['total']/h['total'][0]-1)))
    finite_widths = [w for w in widths if np.isfinite(w)]
    width_min = min(finite_widths) if finite_widths else None
    far = float(h['axion_far'][-1]/target) if h['fv_volume'][-1] > 0 else None
    gates = dict(initial_transfer=bool(np.array_equal(initial, s)),
                 completed=error is None and bool(np.isclose(h['u'][-1], 80.)),
                 energy=drift < 5e-4,
                 wall=bool(widths) and len(finite_widths) == len(widths) and width_min >= 4)
    report = dict(gates=gates, numerical_pass=all(gates.values()), solver_error=error,
                  final_time=float(h['u'][-1]), final_radius=float(h['fv_radius'][-1]),
                  min_radius=float(h['fv_radius'].min()), energy_drift=drift,
                  minimum_wall_cells=width_min, final_far_storage_fraction=far,
                  screening_candidate=all(gates.values()) and h['fv_radius'].min() >= 16 and far is not None and far < .05,
                  interpretation='Pilot only; box isolation, spatial/time convergence and long survival not established.')
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
