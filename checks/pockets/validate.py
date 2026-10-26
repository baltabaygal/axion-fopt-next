"""Predeclared correctness gates for the compiled dilaton+axion pocket model."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from model import (EXECUTABLE, config, run, read_state, read_slice, history, reference, potential,
                   force, dilaton_potential, lambda_for_width, tension, energies, gas_modes, gas_field, grid,
                   radial_modes, bessel_zeros, eigen_state)
from monitor_reference import monitor

GATES = dict(force_difference=2e-7, vacuum_mass_error=2e-6, kink_tension_relative=1e-8,
             kink_residual=1e-12, gas_direct_sum=1e-10, final_ode_error=2e-4, refinement_ratio=3.5,
             snapshot_energy_relative=1e-12, energy_drift_relative=1e-5, bessel_limit_relative=2e-3)
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
    # Radial pocket modes approach the reflective-sphere spectrum, arXiv:2609.24714 Eq. (3.4).
    r0 = 10.;modes = radial_modes(r0, 1., 400., .1, 1.4, 12., h=.01, profile=lambda r: (r > r0).astype(float))
    limit = [abs(m['omega']*r0/bessel_zeros(m['l'], 3)[m['n']]-1) for m in modes if m['l'] <= 2 and m['n'] <= 2]
    details['bessel_limit'] = dict(cases=len(limit), max_relative=float(max(limit)))
    checks['bessel_limit'] = len(limit) == 9 and max(limit) < GATES['bessel_limit_relative']
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
    # Far-field axion energy and all other history columns against the Python measure.
    n, box = 16, 16.;pocket = radial_modes(4., 1., 2., .5, 1.5, 7.)
    s = eigen_state(n, box, 4., 1., pocket, .3, 1, 7.)
    cfg = config(n, box, .025, 2., lam=8., dv=.05, mass_tv=2., every=1, far_margin=1.)
    directory = args.output/'far_field';cfg = run(s, cfg, directory, args.executable);h = history(directory)
    _, final = read_state(directory/'native/final.bin');far_errors = []
    for row, state in ((0, s), (-1, final)):
        e = energies(state, box, 8., .05, 4., far_margin=1.)
        far_errors.append(max(abs(e[key]-h[key][row]) for key in e if key in h.dtype.names)/abs(e['total']))
    details['far_field'] = dict(relative_error=float(max(far_errors)), final_far=float(h['axion_far'][-1]),
                                modes=len(pocket))
    checks['far_field'] = max(far_errors) < GATES['snapshot_energy_relative'] and h['axion_far'][-1] > 0
    # Optional phi^4 coupling: independent force, finite-start ODE and energy
    # comparisons, with timestep refinement separate from spatial effects.
    numerical4 = ((potential(phi+eps, a, LAM, DV, MASS**2, 4)-potential(phi-eps, a, LAM, DV, MASS**2, 4))/(2*eps),
                  (potential(phi, a+eps, LAM, DV, MASS**2, 4)-potential(phi, a-eps, LAM, DV, MASS**2, 4))/(2*eps))
    err4 = float(max(np.max(np.abs(n-f)) for n,f in zip(numerical4, force(phi,a,LAM,DV,MASS**2,4))))
    checks['quartic_force'] = err4 < GATES['force_difference'];details['quartic_force_difference'] = err4
    curvature4=(force(eps, 3., LAM, DV, MASS**2, 4)[0]-force(-eps, 3., LAM, DV, MASS**2, 4)[0])/(2*eps)
    checks['quartic_false_curvature'] = abs(curvature4-LAM/2) < GATES['vacuum_mass_error']
    checks['quartic_true_axion_mass'] = abs(force(1.,.3,LAM,DV,MASS**2,4)[1]/.3-MASS**2) < 1e-12
    errors4=[];s4=small_state(6,6.)
    for dt in (.0025,.00125):
        cfg4=config(6,6.,dt,1.,lam=LAM,dv=DV,mass_tv=MASS,every=1,mass_power=4)
        directory=args.output/f'quartic_{dt}';run(s4,cfg4,directory,args.executable)
        h4=history(directory);_,last4=read_state(directory/'native/final.bin')
        _,first4=read_state(directory/'native/initial.bin')
        err=float(np.max(np.abs(last4-reference(s4,cfg4,[0.,1.])[1])));errors4.append(err)
        checks[f'quartic_{dt}_transfer']=np.array_equal(s4,first4)
        e4=energies(last4,6.,LAM,DV,MASS**2,mass_power=4)
        energy_err=max(abs(e4[k]-h4[k][-1]) for k in e4 if k in h4.dtype.names)/abs(e4['total'])
        drift=float(np.max(np.abs(h4['total']/h4['total'][0]-1)))
        checks[f'quartic_{dt}_energy']=energy_err < GATES['snapshot_energy_relative']
        checks[f'quartic_{dt}_conservation']=drift < GATES['energy_drift_relative']
        # Independent finite-difference Hessian, including mixed curvature,
        # verifies the compiled guard is a coupled-field frequency bound.
        fpp=force(last4[...,0]+eps,last4[...,1],LAM,DV,MASS**2,4)
        fpm=force(last4[...,0]-eps,last4[...,1],LAM,DV,MASS**2,4)
        fap=force(last4[...,0],last4[...,1]+eps,LAM,DV,MASS**2,4)
        fam=force(last4[...,0],last4[...,1]-eps,LAM,DV,MASS**2,4)
        matrix=np.zeros(last4.shape[:-1]+(2,2))
        matrix[...,0,0]=(fpp[0]-fpm[0])/(2*eps)
        matrix[...,0,1]=matrix[...,1,0]=(fap[0]-fam[0])/(2*eps)
        matrix[...,1,1]=(fap[1]-fam[1])/(2*eps)
        expected_bound=dt*np.sqrt(12.+max(0.,np.linalg.eigvalsh(matrix)[...,-1].max()))
        checks[f'quartic_{dt}_hessian_guard']=abs(expected_bound-h4['dt_frequency_bound'][-1])<2e-8
        details[f'quartic_{dt}']=dict(final_ode_error=err,energy_error=energy_err,energy_drift=drift)
    checks['quartic_ode']=max(errors4) < GATES['final_ode_error']
    checks['quartic_refinement']=errors4[0]/errors4[1] > GATES['refinement_ratio']
    details['quartic_refinement_ratio']=errors4[0]/errors4[1]
    pocket4=radial_modes(4.,1.,2.,.5,1.5,7.,mass_power=4)
    s4=eigen_state(16,16.,4.,1.,pocket4,.3,1,7.)
    cfg4=config(16,16.,.0125,2.,lam=8.,dv=.05,mass_tv=2.,every=1,far_margin=1.,mass_power=4)
    directory=args.output/'quartic_far';run(s4,cfg4,directory,args.executable);h4=history(directory)
    _,last4=read_state(directory/'native/final.bin');far_errors4=[]
    for index,state4 in ((0,s4),(-1,last4)):
        e4=energies(state4,16.,8.,.05,4.,far_margin=1.,mass_power=4)
        far_errors4.append(max(abs(e4[k]-h4[k][index]) for k in e4 if k in h4.dtype.names)/abs(e4['total']))
    checks['quartic_far']=max(far_errors4) < GATES['snapshot_energy_relative'] and h4['axion_far'][-1] > 0
    details['quartic_far_error']=max(far_errors4)
    # Fixed sphere, symmetric bond energy and exact semi-discrete current.
    for power in (2,4):
        residuals=[]
        for dt in (.005,.0025):
            s_monitor=small_state(6,6.)
            cfg_monitor=config(6,6.,dt,1.,lam=LAM,dv=DV,mass_tv=MASS,every=1,
                               mass_power=power,monitor_radius=2.2)
            directory=args.output/f'monitor_p{power}_{dt}';run(s_monitor,cfg_monitor,directory,args.executable)
            mh=history(directory);_,last_monitor=read_state(directory/'native/final.bin')
            scale=abs(mh['total'][0]);differences=[]
            for row,state_monitor in ((0,s_monitor),(-1,last_monitor)):
                values=monitor(state_monitor,6.,2.2,LAM,DV,MASS**2,power)
                differences.extend(abs(value-mh[key][row])/scale for key,value in values.items())
            checks[f'monitor_p{power}_{dt}_accounting']=max(differences)<GATES['snapshot_energy_relative']
            integrated=lambda values:np.trapezoid(values,mh['u'])
            residual=np.array([
                mh['monitor_energy'][-1]-mh['monitor_energy'][0]+integrated(mh['monitor_phi_flux']+mh['monitor_axion_flux']),
                mh['monitor_phi_energy'][-1]-mh['monitor_phi_energy'][0]+integrated(mh['monitor_phi_flux']+mh['monitor_exchange']),
                mh['monitor_axion_energy'][-1]-mh['monitor_axion_energy'][0]+integrated(mh['monitor_axion_flux']-mh['monitor_exchange'])])/scale
            residuals.append(float(np.max(np.abs(residual))))
            details[f'monitor_p{power}_{dt}']=dict(snapshot_error=max(differences),continuity_residuals=residual.tolist())
        checks[f'monitor_p{power}_continuity']=max(residuals)<1e-4
        checks[f'monitor_p{power}_refinement']=residuals[0]/residuals[1]>3.5
        details[f'monitor_p{power}_refinement']=residuals[0]/residuals[1]
    # Bad input must fail without creating native output.
    bad = args.output/'bad.in';base = (args.output/'small_0.005/input.in').read_text()
    for label, replacement in [('nonfinite', 'lambda = nan'), ('fractional_grid', 'N = 6.5'),
                               ('expansion', 'expansion = true'), ('negative_mass', 'mass_tv = -1'),
                               ('unsupported_power','mass_power = 3'), ('fractional_power','mass_power = 2.5'),
                               ('negative_monitor','monitor_radius = -1'),('large_monitor','monitor_radius = 3'),
                               ('nonfinite_monitor','monitor_radius = nan')]:
        key = replacement.split(' = ')[0]
        text = '\n'.join(replacement if line.startswith(key+' = ') else line for line in base.splitlines())
        if not any(line.startswith(key+' = ') for line in base.splitlines()):
            text += '\n'+replacement
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
