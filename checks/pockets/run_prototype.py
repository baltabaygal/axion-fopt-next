"""Run the predeclared single-pocket prototype (EXPERIMENT.md) and report gates and outcomes."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from scipy.signal import find_peaks
from model import (ROOT, EXECUTABLE, config, run, history, read_state, read_slice, lambda_for_width,
                   tension, gas_modes, gas_amplitude, pocket_state, radial_collapse, thin_wall, brackets)

WIDTH, DV, R0 = 2., .01, 20.
LAM = lambda_for_width(WIDTH);SIGMA = tension(LAM)
RHO0 = 3*(DV+2*SIGMA/(.85*R0))*.85**4
PARAMETERS = dict(n=128, box=64., dt=.05, r0=R0, wall_width=WIDTH, lambda_=LAM, sigma=SIGMA, delta_v=DV,
                  mass_tv=2., k_min=.15, k_max=.6, seed=0, window_radius=17., window_width=2., rho0=RHO0)
GATES = dict(energy_drift=1e-4, radial_agreement=.25, radial_minimum_radius=8., minimum_wall_cells=4.,
             timestep_radius=.05, spatial_radius=.5, spatial_mean_radius=.25)
OUTCOMES = dict(survival_radius=.4*R0, late_mean_radius=.5*R0, late_window=(200., 300.), leak_fraction=.01,
                control_collapse_volume=.01, control_collapse_time=40., leaky_tv_fraction=.5)
SUITES = dict(control=dict(gas=False, duration=60.), pocket=dict(duration=300., frames=2000),
              timestep=dict(dt=.025, duration=100.), spatial=dict(n=96, duration=100.),
              leaky=dict(mass_tv=.3, duration=100.))


def initial(suite, amplitude_128):
    p = dict(PARAMETERS, **{k: v for k, v in SUITES[suite].items() if k in PARAMETERS})
    if not SUITES[suite].get('gas', True):
        return pocket_state(p['n'], p['box'], R0, WIDTH), 0.
    modes = gas_modes(p['box'], p['k_min'], p['k_max'], p['seed'])
    target = RHO0*4*np.pi/3*R0**3
    # The spatial suite reuses the N=128 amplitude: identical continuum initial data.
    amplitude = amplitude_128 if suite == 'spatial' else gas_amplitude(
        target, p['n'], p['box'], R0, WIDTH, modes, p['window_radius'], p['window_width'], p['mass_tv']**2, LAM, DV)
    return pocket_state(p['n'], p['box'], R0, WIDTH, modes, amplitude, p['window_radius'], p['window_width']), amplitude


def crossings(row, level, dx):
    below = row < level
    edges = np.nonzero(below != np.roll(below, -1))[0]
    nxt = (edges+1) % len(row)
    return np.sort((edges+.5+(level-row[edges])/(row[nxt]-row[edges]))*dx)


def wall_cells(slice_state, dx):
    """10-90% wall thickness in cells along the central slice row; nan unless both walls are seen."""
    row = slice_state[:, len(slice_state)//2, 0]
    inner, outer = crossings(row, .1, dx), crossings(row, .9, dx)
    if len(inner) != 2 or len(outer) != 2 or not outer[0] < inner[0] < inner[1] < outer[1]:
        return np.nan
    return min(inner[0]-outer[0], outer[1]-inner[1])/dx


def gas_diagnostics(state, box):
    """Effective wavenumber and frequency of the axion field inside the false vacuum."""
    dx = box/len(state);phi, a, pa = state[..., 0], state[..., 1], state[..., 3];fv = phi < .5
    grad = sum((np.roll(a, -1, ax)-a)**2 for ax in range(3))/dx**2
    return dict(k_eff=float(np.sqrt(grad[fv].sum()/(a[fv]**2).sum())),
                omega_eff=float(np.sqrt((pa[fv]**2).sum()/(a[fv]**2).sum())),
                dilaton_tv_energy=float(((.5*state[..., 2]**2)[phi > .99]).sum()*dx**3))


def summarize(directory, suite, cfg):
    h = history(directory);native = Path(directory)/'native';dx = cfg['lSide']/cfg['N']
    axion0 = h['axion_kinetic'][0]+h['axion_gradient'][0]+h['coupling'][0]
    walls = []
    for path in sorted(native.glob('slice_*.bin')):
        meta, s = read_slice(path);walls.append((meta['u'], wall_cells(s, dx)))
    walls = np.array(walls)
    drift = float(np.max(np.abs(h['total']-h['total'][0]))/abs(h['total'][0]))
    cov = np.stack([[h['sxx'], h['sxy'], h['sxz']], [h['sxy'], h['syy'], h['syz']], [h['sxz'], h['syz'], h['szz']]])
    eig = np.linalg.eigvalsh(np.moveaxis(cov, -1, 0))
    asphericity = np.where(h['fv_volume'] > 0, (eig[:, 2]-eig[:, 0])/np.maximum(eig.mean(axis=1), 1e-300), np.nan)
    shift = np.sqrt((h['cx']-h['cx'][0])**2+(h['cy']-h['cy'][0])**2+(h['cz']-h['cz'][0])**2)
    alive = h['fv_radius'] >= GATES['radial_minimum_radius']
    report = dict(suite=suite, parameters=cfg, energy_drift=drift, max_dt_frequency=float(h['dt_frequency_bound'].max()),
                  initial_total=float(h['total'][0]), initial_axion=float(axion0),
                  min_radius=float(h['fv_radius'].min()), final_radius=float(h['fv_radius'][-1]),
                  min_wall_cells=float(np.nanmin(walls[walls[:, 0] <= h['u'][alive].max(), 1])) if alive.any() and len(walls) else None,
                  max_tv_fraction=float(np.max(h['axion_tv'])/axion0) if axion0 else None,
                  final_tv_fraction=float(h['axion_tv'][-1]/axion0) if axion0 else None,
                  max_asphericity=float(np.nanmax(asphericity[alive])) if alive.any() else None,
                  max_centroid_shift=float(np.max(shift[h['fv_volume'] > 0])),
                  max_phi2_fv=float(np.max(h['phi2_fv'])))
    report['gates'] = dict(energy_drift=drift < GATES['energy_drift'], frequency_guard=report['max_dt_frequency'] < .5)
    if report['min_wall_cells'] is not None:
        report['gates']['wall_resolution'] = report['min_wall_cells'] >= GATES['minimum_wall_cells']
    return h, report


def compare(h, reference, gate, window=(0., 100.)):
    shared, i, j = np.intersect1d(np.round(h['u'], 9), np.round(reference['u'], 9), return_indices=True)
    keep = (shared >= window[0]) & (shared <= window[1])
    return float(np.max(np.abs(h['fv_radius'][i[keep]]-reference['fv_radius'][j[keep]]))), gate


def main():
    p = argparse.ArgumentParser();p.add_argument('--output', type=Path, required=True)
    p.add_argument('--executable', type=Path, default=EXECUTABLE)
    p.add_argument('--suites', nargs='+', default=list(SUITES), choices=list(SUITES))
    args = p.parse_args();args.output.mkdir(parents=True, exist_ok=False)
    sources = {name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
               for name in ('model.py', 'run_prototype.py', 'pocket_model.h', 'evolve.cpp')}
    contract = dict(parameters=PARAMETERS, gates=GATES, outcomes=OUTCOMES, suites=SUITES, sources=sources,
                    executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                    numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform())
    (args.output/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    _, amplitude_128 = initial('pocket', None)
    reports, histories = {}, {}
    for suite in args.suites:
        spec = dict(PARAMETERS, **SUITES[suite]);started = time.perf_counter()
        state, amplitude = initial(suite, amplitude_128)
        every = int(round(.5/spec['dt']))
        cfg = config(spec['n'], spec['box'], spec['dt'], spec['duration'], lam=LAM, dv=DV, mass_tv=spec['mass_tv'],
                     every=every, frames=spec.get('frames', 0), slices=int(round(2./spec['dt'])))
        try:
            cfg = run(state, cfg, args.output/suite, args.executable);error = None
        except RuntimeError as failure:
            # A solver guard stop still leaves a valid partial history to report.
            cfg = json.loads((args.output/suite/'input.json').read_text());error = str(failure).splitlines()[-1]
        _, transferred = read_state(args.output/suite/'native/initial.bin')
        h, report = summarize(args.output/suite, suite, cfg)
        report['solver_error'] = error;report['completed_until'] = float(h['u'][-1])
        report['gates']['solver_completed'] = error is None
        report['gates']['initial_transfer'] = bool(np.array_equal(transferred, state))
        report['gas_amplitude'] = amplitude;report['elapsed_seconds'] = time.perf_counter()-started
        reports[suite], histories[suite] = report, h
        print(suite, 'radius', round(report['min_radius'], 3), '->', round(report['final_radius'], 3),
              'drift', f"{report['energy_drift']:.2e}", 'seconds', round(report['elapsed_seconds'], 1), flush=True)
    if 'control' in reports:
        h = histories['control'];r = reports['control']
        times, radial = radial_collapse(R0, WIDTH, LAM, DV, 32., .05, 60., .01, 50)
        keep = radial >= GATES['radial_minimum_radius']
        shared, i, j = np.intersect1d(np.round(h['u'], 9), np.round(times[keep], 9), return_indices=True)
        r['radial_difference'] = float(np.max(np.abs(h['fv_radius'][i]-radial[keep][j])))
        r['gates']['radial_reference'] = r['radial_difference'] < GATES['radial_agreement']
        collapsed = h['u'][h['fv_volume'] < OUTCOMES['control_collapse_volume']*h['fv_volume'][0]]
        r['collapse_time'] = float(collapsed[0]) if len(collapsed) else None
        tw_t, tw_r, _ = thin_wall(R0, SIGMA, DV, 0., 60.)
        r['thin_wall_collapse_time'] = float(tw_t[-1])
        resolved = h['fv_radius'] >= GATES['radial_minimum_radius']
        r['energy_drift_resolved'] = float(np.max(np.abs(h['total'][resolved]-h['total'][0]))/abs(h['total'][0]))
        r['outcome'] = dict(collapsed=r['collapse_time'] is not None and r['collapse_time'] <= OUTCOMES['control_collapse_time'])
        np.savetxt(args.output/'control/radial_reference.csv', np.column_stack((times, radial)), delimiter=',',
                   header='u,radius', comments='')
    if 'pocket' in reports:
        h = histories['pocket'];r = reports['pocket'];lo, hi = OUTCOMES['late_window']
        late = (h['u'] >= lo) & (h['u'] <= hi)
        rho0 = r['initial_axion']/(4*np.pi/3*R0**3)
        r['late_mean_radius'] = float(h['fv_radius'][late].mean())
        r['late_radius_std'] = float(h['fv_radius'][late].std())
        r['brackets'] = brackets(R0, SIGMA, DV, rho0, r['initial_total'])
        minima, _ = find_peaks(-h['fv_radius'], prominence=.2)
        r['minima_times'] = h['u'][minima].tolist();r['minima_radii'] = h['fv_radius'][minima].tolist()
        r['breathing_period'] = float(np.diff(h['u'][minima]).mean()) if len(minima) > 1 else None
        tw_t, tw_r, tw_v = thin_wall(R0, SIGMA, DV, rho0, 300.)
        tw_min, _ = find_peaks(-tw_r)
        r['thin_wall'] = dict(first_minimum_time=float(tw_t[tw_min[0]]), first_minimum_radius=float(tw_r[tw_min[0]]),
                              period=float(np.diff(tw_t[tw_min]).mean()), max_speed=float(np.abs(tw_v).max()))
        np.savetxt(args.output/'pocket/thin_wall.csv', np.column_stack((tw_t, tw_r)), delimiter=',',
                   header='u,radius', comments='')
        r['gas'] = {}
        for path in sorted((args.output/'pocket/native').glob('frame_*.bin'))+sorted((args.output/'pocket/native').glob('final.bin')):
            meta, s = read_state(path);r['gas'][f"{meta['u']:g}"] = gas_diagnostics(s, meta['box'])
        r['outcome'] = dict(survives=r['min_radius'] > OUTCOMES['survival_radius'],
                            persists=r['late_mean_radius'] > OUTCOMES['late_mean_radius'],
                            confined=r['final_tv_fraction'] < OUTCOMES['leak_fraction'])
        r['relic_pocket'] = all(r['outcome'].values())
        for suite, gate, name in (('timestep', 'timestep_radius', 'timestep'), ('spatial', 'spatial_radius', 'spatial')):
            if suite in reports:
                diff, tol = compare(histories[suite], h, GATES[gate])
                reports[suite]['radius_difference'] = diff
                reports[suite]['gates'][name+'_agreement'] = diff < tol
        if 'spatial' in reports:
            window = lambda x: x['fv_radius'][(x['u'] >= 50) & (x['u'] <= 100)].mean()
            d = abs(window(histories['spatial'])-window(h))
            reports['spatial']['mean_radius_difference'] = float(d)
            reports['spatial']['gates']['spatial_mean'] = d < GATES['spatial_mean_radius']
        if 'leaky' in reports:
            leaky = histories['leaky'];at = np.argmin(np.abs(leaky['u']-100.));ref = np.argmin(np.abs(h['u']-100.))
            reports['leaky']['outcome'] = dict(
                escaped=reports['leaky']['final_tv_fraction'] > OUTCOMES['leaky_tv_fraction'],
                smaller_than_pocket=bool(leaky['fv_radius'][at] < h['fv_radius'][ref]))
    for suite, report in reports.items():
        report['gates'] = {k: bool(v) for k, v in report['gates'].items()}
        (args.output/suite/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    passed = all(all(r['gates'].values()) for r in reports.values())
    final = dict(numerical_gates_passed=passed, gates={s: r['gates'] for s, r in reports.items()},
                 outcomes={s: r.get('outcome') for s, r in reports.items() if 'outcome' in r},
                 relic_pocket=reports.get('pocket', {}).get('relic_pocket'))
    (args.output/'report.json').write_text(json.dumps(final, indent=2)+'\n')
    print(json.dumps(final, indent=2))
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
