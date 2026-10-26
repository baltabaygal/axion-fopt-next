"""Compiled homogeneous PQ versus the existing fixed-radius cosine reference."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

from model import ROOT, EXECUTABLE, config, run, read_state
sys.path.insert(0, str(ROOT / 'checks'))
from homogeneous import bessel, ode

# Same background, constant mass and finite-start initial data as checks/homogeneous.py.
DURATION = 8.
AMPLITUDES = (.01, .8)
LAMBDAS = (2., 32., 512.)
GATES = dict(ode_field_normalized=2e-4, ode_momentum_normalized=2e-3,
             local_spread_normalized=1e-10, grid_difference_normalized=1e-10,
             refinement_ratio=3., heavy_phase_normalized=2e-3,
             heavy_radial_excursion=2e-3, heavy_improvement_factor=10.)


def polar_reference(times, amplitude, lam):
    theta, velocity = bessel(0., amplitude)
    def rhs(u, z):
        r, v, angle, w = z
        a = 1+u
        return [v, -2*v/a+r*w*w-a*a*(lam*(r*r-1)*r+2*r*(1-np.cos(angle))),
                w, -(2/a+2*v/r)*w-a*a*np.sin(angle)]
    sol = solve_ivp(rhs, (0., DURATION), [1., 0., theta, velocity],
                    method='DOP853', rtol=2e-12, atol=2e-14, dense_output=True,
                    max_step=.01)
    if not sol.success:
        raise RuntimeError(sol.message)
    r, v, angle, w = sol.sol(times)
    a = 1+times
    return np.array([r*np.cos(angle), r*np.sin(angle),
                     a*a*(v*np.cos(angle)-r*w*np.sin(angle)),
                     a*a*(v*np.sin(angle)+r*w*np.cos(angle))]).T


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--executable', type=Path, default=EXECUTABLE)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    contract = dict(gates=GATES, amplitudes=AMPLITUDES, lambdas=LAMBDAS,
                    timesteps=[.00005, .000025], duration=DURATION, N=4, box=2*np.pi,
                    background='tau=1+u, a=tau, H=1/tau^2; conformal time, comoving box',
                    potential='lambda/4*(r^2-1)^2 + r^2*(1-cos(theta)); f_a=m_a=1',
                    initial='r=1, r_prime=0; theta and theta_prime from existing Bessel finite-start data',
                    normalization='errors divided by fixed early-time amplitude; pi_X,Y=a^2 X_prime,Y_prime',
                    stopping='fixed u=8; fail on solver error or any numerical gate',
                    executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                    source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
                                   for f in [Path(__file__), ROOT/'checks/homogeneous.py',
                                             ROOT/'checks/pq_strings/model.py',
                                             ROOT/'checks/pq_strings/pq_model.h',
                                             ROOT/'checks/pq_strings/evolve.cpp']})
    (args.output/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    checks, records, histories = {}, [], {}
    for amplitude in AMPLITUDES:
        phase_errors = []
        for lam in LAMBDAS:
            errors = []
            for dt in (.00005, .000025):
                name = f'A{amplitude:g}_lambda{lam:g}_dt{dt:g}'
                theta0, w0 = bessel(0., amplitude)
                initial = np.zeros((4,4,4,4))
                initial[...] = [np.cos(theta0), np.sin(theta0),
                                -w0*np.sin(theta0), w0*np.cos(theta0)]
                cfg = config(4, 2*np.pi, dt, DURATION, lam=lam, origin=1.,
                             star=1., exponent=0., chi_star=1., every=round(.02/dt))
                directory = args.output/name
                run(initial, cfg, directory, args.executable)
                h = np.genfromtxt(directory/'native/history.csv', delimiter=',', names=True)
                t = h['u']; ref = polar_reference(t, amplitude, lam)
                angle = np.unwrap(np.arctan2(h['y_mean'], h['x_mean']))
                radius = np.hypot(h['x_mean'], h['y_mean'])
                fixed = ode(t, amplitude)[0]
                field_error = float(np.max(np.abs(np.column_stack((h['x_mean'],h['y_mean']))-ref[:,:2]))/amplitude)
                _, final = read_state(directory/'native/final.bin')
                _, transferred = read_state(directory/'native/initial.bin')
                momentum_error = float(np.max(np.abs(final[...,2:]-ref[-1,2:]))/amplitude)
                spread = max(float(np.max(np.ptp(s.reshape(-1,4),axis=0)))/amplitude for s in (transferred, final))
                checks[name+'/transfer'] = bool(np.array_equal(initial,transferred))
                checks[name+'/ode_field'] = field_error < GATES['ode_field_normalized']
                checks[name+'/ode_momentum'] = momentum_error < GATES['ode_momentum_normalized']
                checks[name+'/spread'] = spread < GATES['local_spread_normalized']
                checks[name+'/gradient'] = bool(np.max(np.abs(h['gradient'])) < 1e-25)
                checks[name+'/time'] = bool(np.max(abs(h['tau']-(1+t))) < 1e-12 and np.max(abs(h['a']-(1+t))) < 1e-12 and np.all(h['chi']==1))
                rec = dict(name=name, amplitude=amplitude, lam=lam, radial_mass=float(np.sqrt(2*lam)), dt=dt,
                           ode_field_error=field_error, ode_momentum_error=momentum_error,
                           phase_difference_normalized=float(np.max(abs(angle-fixed))/amplitude),
                           max_radial_excursion=float(np.max(abs(radius-1))), local_spread_normalized=spread)
                records.append(rec); errors.append(field_error); histories[name] = h
                np.savetxt(directory/'comparison.csv',np.column_stack((t,angle,fixed,radius,ref)),delimiter=',',
                           header='u,theta_pq,theta_fixed_radius,r_pq,X_reference,Y_reference,pi_X_reference,pi_Y_reference',comments='')
            ratio = errors[0]/errors[1]
            checks[f'A{amplitude}/lambda{lam}/refinement'] = ratio > GATES['refinement_ratio']
            records[-1]['refinement_ratio'] = ratio
            phase_errors.append(records[-1]['phase_difference_normalized'])
        checks[f'A{amplitude}/heavy_phase'] = phase_errors[-1] < GATES['heavy_phase_normalized']
        checks[f'A{amplitude}/heavy_radius'] = records[-1]['max_radial_excursion'] < GATES['heavy_radial_excursion']
        checks[f'A{amplitude}/mass_trend'] = bool(np.all(np.diff(phase_errors)<0))
        checks[f'A{amplitude}/heavy_improvement'] = phase_errors[0]/phase_errors[-1] > GATES['heavy_improvement_factor']
    # A larger grid with identical homogeneous data is an invariance check, not spatial convergence.
    base = args.output/records[-1]['name']
    cfg = json.loads((base/'input.json').read_text()); cfg['N']=8
    _, initial = read_state(base/'native/initial.bin')
    run(np.broadcast_to(initial[0,0,0], (8,8,8,4)).copy(),cfg,args.output/'grid_N8',args.executable)
    grid = np.genfromtxt(args.output/'grid_N8/native/history.csv',delimiter=',',names=True)
    h = histories[records[-1]['name']]
    checks['grid_invariance'] = bool(max(np.max(abs(grid[k]-h[k])) for k in ('x_mean','y_mean','r_mean'))/AMPLITUDES[-1] < GATES['grid_difference_normalized'])
    report = dict(passed=all(checks.values()),checks={k:bool(v) for k,v in checks.items()},runs=records)
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
