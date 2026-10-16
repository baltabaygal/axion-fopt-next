"""Independent references, initial states and binary I/O for the pocket model.

Static background, dilaton vev v=1 and unit axion decay constant. Fields are
phi (dilaton: false vacuum 0, true vacuum 1) and a (axion, m_a^2=mass_tv^2*phi^2).
"""
from pathlib import Path
import json
import struct
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'build/pockets/pocket_evolve'


def dilaton_potential(phi, lam, dv):
    """Expanded form; V(1)=0 and V(0)=dv. The compiled model uses a factored form."""
    return lam/4*(phi**2-2*phi**3+phi**4)+dv*(1-4*phi**3+3*phi**4)


def potential(phi, a, lam, dv, m2):
    return dilaton_potential(phi, lam, dv)+.5*m2*phi**2*a**2


def force(phi, a, lam, dv, m2):
    """(dV/dphi, dV/da), differentiated from the expanded polynomial."""
    return (lam/4*(2*phi-6*phi**2+4*phi**3)+dv*(12*phi**3-12*phi**2)+m2*phi*a**2,
            m2*phi**2*a)


def lambda_for_width(width):
    """Thin-wall kink phi=(1+tanh(z/width))/2 requires lambda=8/width^2."""
    return 8/width**2


def tension(lam):
    return np.sqrt(lam/2)/6


def grid(n, box):
    return (np.arange(n)+.5)*box/n


def radius(n, box, centre=None):
    """Periodic nearest-image distance of every cell centre from centre."""
    centre = np.full(3, box/2) if centre is None else np.asarray(centre, float)
    x = grid(n, box)
    d = [np.abs(x-c) for c in centre]
    d = [np.minimum(v, box-v) for v in d]
    return np.sqrt(d[0][:, None, None]**2+d[1][None, :, None]**2+d[2][None, None, :]**2)


def energies(state, box, lam, dv, m2):
    """Volume-integrated components and region split, matching the driver's definitions."""
    n = len(state);dx = box/n;cell = dx**3
    phi, a, pphi, pa = (state[..., i] for i in range(4))
    gphi = sum((np.roll(phi, -1, ax)-phi)**2 for ax in range(3))/(2*dx*dx)
    ga = sum((np.roll(a, -1, ax)-a)**2 for ax in range(3))/(2*dx*dx)
    coupling = .5*m2*phi**2*a**2
    axion = .5*pa**2+ga+coupling
    fv = phi < .5;tv = phi > .99
    out = dict(dilaton_kinetic=.5*np.sum(pphi**2), dilaton_gradient=np.sum(gphi),
               dilaton_potential=np.sum(dilaton_potential(phi, lam, dv)),
               axion_kinetic=.5*np.sum(pa**2), axion_gradient=np.sum(ga), coupling=np.sum(coupling),
               fv_volume=float(np.count_nonzero(fv)), axion_fv=np.sum(axion[fv]),
               axion_wall=np.sum(axion[~fv & ~tv]), axion_tv=np.sum(axion[tv]))
    out = {k: v*cell for k, v in out.items()}
    out['total'] = sum(out[k] for k in ('dilaton_kinetic', 'dilaton_gradient', 'dilaton_potential',
                                          'axion_kinetic', 'axion_gradient', 'coupling'))
    out['fv_radius'] = np.cbrt(3*out['fv_volume']/(4*np.pi))
    return out


def dilaton_sphere(n, box, r0, width):
    """False-vacuum ball of radius r0 in true vacuum, wall at rest."""
    return .5*(1+np.tanh((radius(n, box)-r0)/width))


def gas_modes(box, k_min, k_max, seed):
    """Fourier modes 2*pi*m/box in a shell, with unit-temperature equipartition amplitudes.

    Field amplitudes have variance 1/omega^2 and momenta variance 1, independent,
    so each mode carries equal mean energy and an equal mix of propagation senses.
    """
    m_max = int(np.floor(k_max*box/(2*np.pi)))
    m = np.arange(-m_max, m_max+1)
    vectors = np.stack(np.meshgrid(m, m, m, indexing='ij'), -1).reshape(-1, 3)
    k = 2*np.pi*np.linalg.norm(vectors, axis=1)/box
    keep = (k >= k_min) & (k <= k_max)
    vectors, k = vectors[keep], k[keep]
    rng = np.random.default_rng(seed)
    field = (rng.normal(size=len(k))+1j*rng.normal(size=len(k)))/k
    momentum = rng.normal(size=len(k))+1j*rng.normal(size=len(k))
    return dict(vectors=vectors, k=k, field=field, momentum=momentum, box=box)


def gas_field(modes, n):
    """Evaluate Re sum_m c_m exp(i k_m.x) at cell centres; grid-independent continuum function."""
    box = modes['box'];dx = box/n
    if np.abs(modes['vectors']).max() >= n/2:
        raise ValueError('Gas modes exceed the grid Nyquist frequency')
    k = 2*np.pi*modes['vectors']/box
    shift = np.exp(1j*k.sum(axis=1)*dx/2)
    index = tuple((modes['vectors'] % n).T)
    fields = []
    for coefficient in (modes['field'], modes['momentum']):
        spectrum = np.zeros((n, n, n), complex)
        np.add.at(spectrum, index, coefficient*shift)
        fields.append(np.fft.ifftn(spectrum).real*n**3)
    return fields


def window(n, box, r_window, width):
    return .5*(1-np.tanh((radius(n, box)-r_window)/width))


def pocket_state(n, box, r0, wall_width, modes=None, amplitude=0., window_radius=None, window_width=None):
    s = np.zeros((n, n, n, 4))
    s[..., 0] = dilaton_sphere(n, box, r0, wall_width)
    if modes is not None and amplitude:
        a, pa = gas_field(modes, n)
        w = window(n, box, window_radius, window_width)
        s[..., 1] = amplitude*w*a
        s[..., 3] = amplitude*w*pa
    return s


def gas_amplitude(target_energy, n, box, r0, wall_width, modes, window_radius, window_width, m2, lam, dv):
    """Amplitude giving the requested lattice axion energy (all axion terms are quadratic in a)."""
    s = pocket_state(n, box, r0, wall_width, modes, 1., window_radius, window_width)
    e = energies(s, box, lam, dv, m2)
    unit = e['axion_kinetic']+e['axion_gradient']+e['coupling']
    return float(np.sqrt(target_energy/unit))


def laplacian(field, dx):
    return sum(np.roll(field, 1, ax)+np.roll(field, -1, ax)-2*field for ax in range(3))/dx**2


def reference(initial, cfg, times):
    """Independent semi-discrete DOP853 evolution on the same periodic stencil."""
    shape = initial.shape;dx = cfg['lSide']/shape[0]
    lam, dv, m2 = cfg['lambda'], cfg['delta_v'], cfg['mass_tv']**2

    def rhs(_, flat):
        s = flat.reshape(shape)
        fphi, fa = force(s[..., 0], s[..., 1], lam, dv, m2)
        acceleration = np.stack((laplacian(s[..., 0], dx)-fphi, laplacian(s[..., 1], dx)-fa), -1)
        return np.concatenate((s[..., 2:], acceleration), -1).ravel()
    result = solve_ivp(rhs, (0, max(times)), initial.ravel(), t_eval=times, method='DOP853',
                       rtol=2e-10, atol=2e-12, max_step=.02)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y.T.reshape((len(times),)+shape)


def radial_collapse(r0, width, lam, dv, r_max, dr, duration, dt, every):
    """Spherically symmetric dilaton-only reference: phi_tt = r^-2 (r^2 phi_r)_r - V'(phi).

    Cell-centred r grid, regular at the origin, phi=1 held at r_max. Returns
    times and the radius where phi crosses 1/2 (outermost crossing).
    """
    r = (np.arange(int(round(r_max/dr)))+.5)*dr
    faces = np.arange(len(r)+1)*dr
    phi = .5*(1+np.tanh((r-r0)/width));velocity = np.zeros_like(r)

    def acceleration(f):
        padded = np.concatenate((f[:1], f, [1.]))
        flux = faces**2*np.diff(padded)/dr
        return np.diff(flux)/(r**2*dr)-force(f, 0., lam, dv, 0.)[0]

    def crossing(f):
        inside = np.nonzero(f < .5)[0]
        if not len(inside):
            return 0.
        i = inside[-1]
        return r[i]+(.5-f[i])/(f[i+1]-f[i])*dr
    steps = int(round(duration/dt));times, radii = [], []
    accel = acceleration(phi)
    for step in range(steps+1):
        if step % every == 0:
            times.append(step*dt);radii.append(crossing(phi))
        velocity += .5*dt*accel;phi += dt*velocity
        accel = acceleration(phi);velocity += .5*dt*accel
    return np.array(times), np.array(radii)


def thin_wall(r0, sigma, dv, rho0, duration, max_step=.05):
    """Relativistic thin wall with an adiabatic isotropic gas, P=rho0*(r0/R)^4/3.

    Momentum p=sigma*R^2*gamma*Rdot; d p/dt = -2 sigma R/gamma - dv R^2 + P R^2.
    """
    def rhs(_, y):
        r, p = y;m = sigma*r*r;g = np.sqrt(1+(p/m)**2)
        return [p/(m*g), -2*sigma*r/g-dv*r*r+rho0*(r0/r)**4/3*r*r]
    stop = lambda _, y: y[0]-.05*r0
    stop.terminal = True
    result = solve_ivp(rhs, (0, duration), [r0, 0.], max_step=max_step, events=stop, rtol=1e-10, atol=1e-10)
    r, p = result.y;m = sigma*r*r
    return result.t, r, p/(m*np.sqrt(1+(p/m)**2))


def brackets(r0, sigma, dv, rho0, energy):
    """Equilibrium radii: reversible adiabatic gas, and all released energy kept as gas heat.

    The second solves 3*sigma*A + 4*dv*V = energy, using pressure balance
    rho/3 = dv + 2*sigma/R and wall energy sigma*A.
    """
    adiabatic = brentq(lambda r: rho0*(r0/r)**4/3-dv-2*sigma/r, 1e-3*r0, 10*r0)
    dumped = brentq(lambda r: 3*sigma*4*np.pi*r*r+4*dv*4*np.pi/3*r**3-energy, 1e-3*r0, 10*r0)
    return dict(adiabatic=adiabatic, energy_conserving=dumped)


def write_initial(path, state, box):
    state = np.asarray(state, dtype='<f8');n = state.shape[0]
    if state.shape != (n, n, n, 4) or not np.isfinite(state).all() or box <= 0:
        raise ValueError('Require finite cubic (N,N,N,4) state and positive box')
    with Path(path).open('wb') as f:
        f.write(b'PKIN0001'+struct.pack('<Id', n, box))
        f.write(state.tobytes())


def _read(path, magic, cells):
    raw = Path(path).read_bytes()
    if raw[:8] != magic or len(raw) < 28:
        raise ValueError('Invalid magic/header')
    n, box, u = struct.unpack('<I2d', raw[8:28])
    if not 4 <= n <= 128 or not np.isfinite([box, u]).all() or box <= 0:
        raise ValueError('Invalid metadata')
    if len(raw) != 28+32*cells(n):
        raise ValueError('Invalid byte count')
    data = np.frombuffer(raw, dtype='<f8', offset=28).copy()
    if not np.isfinite(data).all():
        raise ValueError('Nonfinite data')
    return dict(n=n, box=box, u=u), data


def read_state(path):
    meta, data = _read(path, b'PKST0001', lambda n: n**3)
    return meta, data.reshape(meta['n'], meta['n'], meta['n'], 4)


def read_slice(path):
    """Plane z index n/2: (n, n, 4) array of phi, a, pi_phi, pi_a."""
    meta, data = _read(path, b'PKSL0001', lambda n: n**2)
    return meta, data.reshape(meta['n'], meta['n'], 4)


def config(n, box, dt, duration, *, lam, dv, mass_tv, every=10, frames=0, slices=0):
    return dict(N=n, lSide=box, dt=dt, tMax=duration, expansion=False, evolver='PV2',
                **{'lambda': lam}, delta_v=dv, mass_tv=mass_tv,
                output_every=every, snapshot_every=frames, slice_every=slices)


def run(initial, cfg, root, executable=EXECUTABLE):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    cfg = dict(cfg)
    write_initial(root/'input.bin', initial, cfg['lSide'])
    cfg.update(initial_state=str((root/'input.bin').resolve()), result_dir=str((root/'native').resolve()))
    (root/'input.json').write_text(json.dumps(cfg, indent=2)+'\n')
    text = '\n'.join(f'{k} = {str(v).lower() if isinstance(v, bool) else v}' for k, v in cfg.items())+'\n'
    (root/'input.in').write_text(text)
    result = subprocess.run([str(executable), 'input='+str((root/'input.in').resolve())],
                            capture_output=True, text=True)
    (root/'run.log').write_text(result.stdout+result.stderr)
    if result.returncode:
        raise RuntimeError(f'Pocket solver failed: {root}/run.log\n{result.stderr[-1500:]}')
    return cfg


def history(directory):
    return np.genfromtxt(Path(directory)/'native/history.csv', delimiter=',', names=True)
