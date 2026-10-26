"""Run the predeclared pocket lifetime and leak-rate experiment (LEAK_EXPERIMENT.md)."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from model import (EXECUTABLE, config, run, history, read_state, read_slice, lambda_for_width, tension,
                   radial_modes, eigen_state, energies, radius, grid)
from run_prototype import wall_cells, gas_diagnostics

WIDTH, DV, R0, BOX = 2., .01, 20., 64.
LAM = lambda_for_width(WIDTH);SIGMA = tension(LAM)
GAS = dict(k_min=.15, k_max=.6, r_max=31., h=.01, seed=0, r_cut=30., energy=.0464*4*np.pi/3*R0**3)
RUNS = {
    'control': dict(n=128, dt=.025, mass_tv=2., duration=30., gas=False),
    'm2_N128': dict(n=128, dt=.025, mass_tv=2., duration=400., frames=100.),
    'm2_N128_dt': dict(n=128, dt=.0125, mass_tv=2., duration=100.),
    'm2_N96': dict(n=96, dt=.025, mass_tv=2., duration=400.),
    'm2_N192': dict(n=192, dt=.025, mass_tv=2., duration=400.),
    'm1.5_N128': dict(n=128, dt=.025, mass_tv=1.5, duration=400.),
    'm3_N128': dict(n=128, dt=.025, mass_tv=3., duration=400.),
    'm4_N128': dict(n=128, dt=.025, mass_tv=4., duration=400.),
}
GATES = dict(initial_far=1e-6, energy_drift=5e-4, survival_radius=.4*R0, timestep_radius=.05,
             timestep_far=.005, minimum_wall_cells=4.)
CLASSES = dict(leak_window=(25., 125.), minimum_window=25., resolution_relative=.2, no_leak=1e-3,
               steep=.1, weak=.5, collapse_volume=.01)


def pocket_modes(mass_tv):
    return radial_modes(R0, WIDTH, mass_tv, GAS['k_min'], GAS['k_max'], GAS['r_max'], GAS['h'])


def amplitude(mass_tv, modes):
    """Fixed on N=128 so every grid samples one continuum initial state."""
    unit = eigen_state(128, BOX, R0, WIDTH, modes, 1., GAS['seed'], GAS['r_cut'])
    e = energies(unit, BOX, LAM, DV, mass_tv**2)
    return float(np.sqrt(GAS['energy']/(e['axion_kinetic']+e['axion_gradient']+e['coupling'])))


def far_frequency(state, box, margin=6.):
    """Effective axion frequency in the far field, to separate propagating escapees from tails."""
    n = len(state);dx = box/n;phi, a, pa = state[..., 0], state[..., 1], state[..., 3];fv = phi < .5
    if not fv.any():
        return None
    x = grid(n, box);centre = [x[i].mean() for i in np.nonzero(fv)]
    far = radius(n, box, centre) > np.cbrt(3*fv.sum()*dx**3/(4*np.pi))+margin
    return float(np.sqrt((pa[far]**2).sum()/(a[far]**2).sum())) if (a[far]**2).sum() > 0 else None


def execute(name, root, executable, modes, amplitudes):
    spec = RUNS[name];started = time.perf_counter()
    gas = spec.get('gas', True);m = spec['mass_tv']
    state = eigen_state(spec['n'], BOX, R0, WIDTH, modes[m] if gas else None, amplitudes[m] if gas else 0.,
                        GAS['seed'], GAS['r_cut'])
    cfg = config(spec['n'], BOX, spec['dt'], spec['duration'], lam=LAM, dv=DV, mass_tv=m,
                 every=int(round(.5/spec['dt'])), slices=int(round(4./spec['dt'])),
                 frames=int(round(spec['frames']/spec['dt'])) if 'frames' in spec else 0)
    try:
        cfg = run(state, cfg, root/name, executable);error = None
    except RuntimeError as failure:
        cfg = json.loads((root/name/'input.json').read_text());error = str(failure).splitlines()[-1]
    _, transferred = read_state(root/name/'native/initial.bin')
    print(name, 'done', round(time.perf_counter()-started, 1), 's', error or '', flush=True)
    return dict(cfg=cfg, error=error, transfer=bool(np.array_equal(transferred, state)),
                amplitude=amplitudes[m] if gas else 0., elapsed=time.perf_counter()-started)


def summarize(name, root, result):
    h = history(root/name);spec = RUNS[name];dx = BOX/spec['n']
    axion0 = h['axion_kinetic'][0]+h['axion_gradient'][0]+h['coupling'][0]
    dead = np.nonzero(h['fv_radius'] < GATES['survival_radius'])[0]
    life = float(h['u'][dead[0]]) if len(dead) else None
    alive = h['u'] < (life if life is not None else np.inf)
    end = float(h['u'][-1])
    r = dict(name=name, spec=spec, parameters=result['cfg'], solver_error=result['error'], completed_until=end,
             elapsed_seconds=result['elapsed'], gas_amplitude=result['amplitude'], lifetime=life,
             censored=life is None, initial_total=float(h['total'][0]), initial_axion=float(axion0),
             energy_drift_alive=float(np.max(np.abs(h['total'][alive]-h['total'][0]))/abs(h['total'][0])),
             max_dt_frequency=float(h['dt_frequency_bound'].max()),
             min_radius=float(h['fv_radius'].min()), final_radius=float(h['fv_radius'][-1]))
    walls = [(meta['u'], wall_cells(s, dx)) for meta, s in map(read_slice, sorted((root/name/'native').glob('slice_*.bin')))]
    alive_walls = [w for u, w in walls if life is None or u < life]
    r['min_wall_cells_alive'] = float(np.nanmin(alive_walls)) if alive_walls and not np.all(np.isnan(alive_walls)) else None
    r['gates'] = dict(initial_transfer=result['transfer'],
                      guard_not_while_alive=result['error'] is None or (life is not None and end >= life))
    if axion0 > 0:
        far, tv = h['axion_far']/axion0, h['axion_tv']/axion0
        r['f_far'] = {f'{t:g}': float(far[np.argmin(np.abs(h['u']-t))]) for t in (0, 25, 50, 100, 200, 300, 400) if t <= end}
        r['f_tv'] = {f'{t:g}': float(tv[np.argmin(np.abs(h['u']-t))]) for t in (0, 25, 50, 100, 200, 300, 400) if t <= end}
        lo, hi = CLASSES['leak_window'];hi = min(hi, end, life if life is not None else np.inf)
        window = (h['u'] >= lo) & (h['u'] <= hi)
        r['leak_window'] = [lo, float(hi)]
        r['leak_rate'] = float(np.polyfit(h['u'][window], far[window], 1)[0]) if hi-lo >= CLASSES['minimum_window'] else None
        tail = alive & (h['u'] >= h['u'][alive].max()-100)
        r['late_dR_dt'] = float(np.polyfit(h['u'][tail], h['fv_radius'][tail], 1)[0]) if tail.sum() > 10 else None
        r['gates'].update(initial_confinement=far[0] < GATES['initial_far'],
                          energy_drift=r['energy_drift_alive'] < GATES['energy_drift'])
        if name == 'm2_N192':
            r['gates']['wall_resolution'] = r['min_wall_cells_alive'] is not None and r['min_wall_cells_alive'] >= GATES['minimum_wall_cells']
        r['gas'] = {}
        for path in sorted((root/name/'native').glob('frame_*.bin'))+sorted((root/name/'native').glob('final.bin')):
            meta, s = read_state(path)
            r['gas'][f"{meta['u']:g}"] = dict(gas_diagnostics(s, BOX), far_omega_eff=far_frequency(s, BOX))
    else:
        collapsed = h['u'][h['fv_volume'] < CLASSES['collapse_volume']*h['fv_volume'][0]]
        r['collapse_time'] = float(collapsed[0]) if len(collapsed) else None
    r['gates'] = {k: bool(v) for k, v in r['gates'].items()}
    (root/name/'summary.json').write_text(json.dumps(r, indent=2)+'\n')
    return r, h


def main():
    p = argparse.ArgumentParser();p.add_argument('--output', type=Path, required=True)
    p.add_argument('--executable', type=Path, default=EXECUTABLE);p.add_argument('--jobs', type=int, default=4)
    args = p.parse_args();args.output.mkdir(parents=True, exist_ok=False)
    sources = {name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
               for name in ('model.py', 'run_leak.py', 'run_prototype.py', 'pocket_model.h', 'evolve.cpp')}
    masses = sorted({spec['mass_tv'] for spec in RUNS.values()})
    modes = {m: pocket_modes(m) for m in masses}
    amplitudes = {m: amplitude(m, modes[m]) for m in masses}
    contract = dict(runs=RUNS, gas=GAS, gates=GATES, classes=CLASSES, wall_width=WIDTH, delta_v=DV, r0=R0, box=BOX,
                    lambda_=LAM, sigma=SIGMA, bound_states={f'{m:g}': sum(2*x['l']+1 for x in modes[m]) for m in masses},
                    amplitudes={f'{m:g}': a for m, a in amplitudes.items()}, sources=sources,
                    executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                    numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform())
    (args.output/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    order = sorted(RUNS, key=lambda k: -RUNS[k]['n']**3*RUNS[k]['duration']/RUNS[k]['dt'])
    with ThreadPoolExecutor(args.jobs) as pool:
        results = dict(zip(order, pool.map(lambda k: execute(k, args.output, args.executable, modes, amplitudes), order)))
    reports, histories = {}, {}
    for name in RUNS:
        reports[name], histories[name] = summarize(name, args.output, results[name])
    tc = reports['control'].get('collapse_time')
    for r in reports.values():
        if 'leak_rate' in r and tc:
            r['lifetime_over_collapse'] = (r['lifetime'] if r['lifetime'] is not None else r['completed_until'])/tc
    a, b = histories['m2_N128'], histories['m2_N128_dt']
    shared, i, j = np.intersect1d(np.round(a['u'], 9), np.round(b['u'], 9), return_indices=True)
    keep = shared <= 100
    dt_radius = float(np.max(np.abs(a['fv_radius'][i[keep]]-b['fv_radius'][j[keep]])))
    dt_far = abs(reports['m2_N128']['f_far']['100']-reports['m2_N128_dt']['f_far']['100'])
    reports['m2_N128_dt']['gates'].update(timestep_radius=dt_radius < GATES['timestep_radius'],
                                          timestep_far=dt_far < GATES['timestep_far'])
    reports['m2_N128_dt'].update(timestep_radius_difference=dt_radius, timestep_far_difference=dt_far)
    fine, base = reports['m2_N192'], reports['m2_N128']
    if fine['f_far'].get('100', 0) < CLASSES['no_leak'] and base['f_far'].get('100', 0) < CLASSES['no_leak']:
        resolution = 'no measurable leak'
    else:
        rel = CLASSES['resolution_relative']
        same_rate = None not in (fine['leak_rate'], base['leak_rate']) and abs(fine['leak_rate']-base['leak_rate']) <= rel*abs(base['leak_rate'])
        same_far = abs(fine['f_far']['100']-base['f_far']['100']) <= rel*base['f_far']['100']
        resolution = 'resolution-independent leak' if same_rate and same_far else 'resolution-dependent'
    low, high = reports['m1.5_N128']['leak_rate'], reports['m4_N128']['leak_rate']
    ratio = high/low if None not in (low, high) and low > 0 else None
    mass = None if ratio is None else ('steep' if ratio <= CLASSES['steep'] else 'weak' if ratio >= CLASSES['weak'] else 'moderate')
    for name, r in reports.items():
        (args.output/name/'summary.json').write_text(json.dumps(r, indent=2)+'\n')
    passed = all(all(r['gates'].values()) for r in reports.values())
    final = dict(numerical_gates_passed=passed, gates={k: r['gates'] for k, r in reports.items()},
                 control_collapse_time=tc,
                 lifetimes={k: dict(lifetime=r['lifetime'], censored=r['censored'], over_collapse=r.get('lifetime_over_collapse'))
                            for k, r in reports.items() if 'leak_rate' in r},
                 leak_rates={k: r['leak_rate'] for k, r in reports.items() if 'leak_rate' in r},
                 f_far_100={k: r['f_far'].get('100') for k, r in reports.items() if 'f_far' in r},
                 resolution_class=resolution, mass_ratio=ratio, mass_class=mass)
    (args.output/'report.json').write_text(json.dumps(final, indent=2)+'\n')
    print(json.dumps(final, indent=2))
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
