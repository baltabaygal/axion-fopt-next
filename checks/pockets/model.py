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
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import sph_harm_y, spherical_jn
import importlib.util
# References are also loaded directly by file path, outside this directory's
# import path (e.g. optional integration tests).
_wall_spec=importlib.util.spec_from_file_location('pocket_structured_wall',Path(__file__).with_name('structured_wall.py'))
structured_wall=importlib.util.module_from_spec(_wall_spec)
_wall_spec.loader.exec_module(structured_wall)

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'build/pockets/pocket_evolve'
MAX_N = 384


def dilaton_potential(phi, lam, dv, vacuum_stiffening=0.):
    """Expanded form; V(1)=0 and V(0)=dv. The compiled model uses a factored form."""
    if vacuum_stiffening:
        return structured_wall.potential(phi,lam,dv,vacuum_stiffening)
    return lam/4*(phi**2-2*phi**3+phi**4)+dv*(1-4*phi**3+3*phi**4)


def _check_power(mass_power):
    if mass_power not in (2, 4):
        raise ValueError('mass_power must be 2 or 4')


def potential(phi, a, lam, dv, m2, mass_power=2, vacuum_stiffening=0.):
    _check_power(mass_power)
    return dilaton_potential(phi, lam, dv,vacuum_stiffening)+.5*m2*phi**mass_power*a**2


def force(phi, a, lam, dv, m2, mass_power=2, vacuum_stiffening=0.):
    """(dV/dphi, dV/da), differentiated from the expanded polynomial."""
    _check_power(mass_power)
    if vacuum_stiffening:
        return (structured_wall.force(phi,lam,dv,vacuum_stiffening)+.5*mass_power*m2*phi**(mass_power-1)*a**2,
                m2*phi**mass_power*a)
    return (lam/4*(2*phi-6*phi**2+4*phi**3)+dv*(12*phi**3-12*phi**2)
            +.5*mass_power*m2*phi**(mass_power-1)*a**2,
            m2*phi**mass_power*a)


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


def energies(state, box, lam, dv, m2, far_margin=6., mass_power=2, vacuum_stiffening=0.):
    """Volume-integrated components and region split, matching the driver's definitions.

    'axion_far' counts axion energy beyond R_eff+far_margin from the false-vacuum
    centroid (periodic nearest image), R_eff the volume-equivalent radius.
    """
    _check_power(mass_power)
    n = len(state);dx = box/n;cell = dx**3
    phi, a, pphi, pa = (state[..., i] for i in range(4))
    gphi = sum((np.roll(phi, -1, ax)-phi)**2 for ax in range(3))/(2*dx*dx)
    ga = sum((np.roll(a, -1, ax)-a)**2 for ax in range(3))/(2*dx*dx)
    coupling = .5*m2*phi**mass_power*a**2
    axion = .5*pa**2+ga+coupling
    fv = phi < .5;tv = phi > .99
    far = np.zeros_like(fv)
    if fv.any():
        x = grid(n, box);idx = np.nonzero(fv)
        centre = [x[i].mean() for i in idx]
        far = radius(n, box, centre) > np.cbrt(3*fv.sum()*cell/(4*np.pi))+far_margin
    out = dict(dilaton_kinetic=.5*np.sum(pphi**2), dilaton_gradient=np.sum(gphi),
               dilaton_potential=np.sum(dilaton_potential(phi, lam, dv,vacuum_stiffening)),
               axion_kinetic=.5*np.sum(pa**2), axion_gradient=np.sum(ga), coupling=np.sum(coupling),
               fv_volume=float(np.count_nonzero(fv)), axion_fv=np.sum(axion[fv]),
               axion_wall=np.sum(axion[~fv & ~tv]), axion_tv=np.sum(axion[tv]), axion_far=np.sum(axion[far]))
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


def radial_modes(r0, width, mass_tv, k_min, k_max, r_max, h=.01, profile=None, mass_power=2):
    """Bound pocket modes with k_min <= omega <= k_max for a static spherical wall.

    Solves -u'' + [l(l+1)/r^2 + mass_tv^2 phi(r)^2] u = omega^2 u, u = r R(r),
    u(0)=u(r_max)=0, by second-order finite differences. The default wall is
    phi=[1+tanh((r-r0)/width)]/2. In the Dirichlet limit (large mass, sharp wall)
    omega_nl -> alpha_nl/r0, arXiv:2609.24714 Eq. (3.4). Each u has int u^2 dr = 1.
    """
    _check_power(mass_power)
    r = h*np.arange(1, int(round(r_max/h)))
    phi = .5*(1+np.tanh((r-r0)/width)) if profile is None else profile(r)
    off = -np.ones(len(r)-1)/h**2;modes = []
    for l in range(10000):
        diagonal = 2/h**2+l*(l+1)/r**2+mass_tv**2*phi**mass_power
        w2, u = eigh_tridiagonal(diagonal, off, select='v', select_range=(k_min**2, k_max**2))
        if not len(w2):
            if l*(l+1) > (k_max*r0)**2:
                break
            continue
        for n, (value, vector) in enumerate(zip(w2, u.T)):
            modes.append(dict(l=l, n=n, omega=float(np.sqrt(value)), r=r, u=vector/np.sqrt(h)))
    return modes


def bessel_zeros(l, count):
    """First positive zeros alpha_l1..alpha_l,count of the spherical Bessel function j_l."""
    x = np.linspace(.5, 3*count+l+10, 20000);f = lambda v: spherical_jn(l, v);y = f(x)
    brackets = np.nonzero(np.sign(y[:-1]) != np.sign(y[1:]))[0][:count]
    return np.array([brentq(f, x[i], x[i+1], xtol=1e-14) for i in brackets])


def real_sph_harm(l, m, theta, azimuth):
    """Orthonormal real spherical harmonics; theta is the polar angle."""
    y = sph_harm_y(l, abs(m), theta, azimuth)
    if m > 0:
        return np.sqrt(2)*(-1)**m*y.real
    if m < 0:
        return np.sqrt(2)*(-1)**m*y.imag
    return y.real


def eigen_gas(modes, n, box, seed, r_cut=None, chunk_size=65536):
    """Unit-temperature equipartition gas of pocket standing waves R_nl(r) Y_lm.

    Field amplitudes have variance 1/omega^2 and momenta variance 1, independent.
    Random numbers are drawn in (l, n, m) order, so a seed fixes the continuum
    field independently of the grid. Cells beyond r_cut are left at zero.
    """
    if chunk_size<1:
        raise ValueError('chunk_size must be positive')
    x = grid(n, box)-box/2;rng = np.random.default_rng(seed)
    groups=[]
    for l in sorted({mode['l'] for mode in modes}):
        entries=[]
        for mode in (mode for mode in modes if mode['l'] == l):
            coefficients=[(rng.normal()/mode['omega'],rng.normal()) for _ in range(2*l+1)]
            entries.append((mode,coefficients))
        groups.append((l,entries))
    # Keep the exact coefficient and per-cell accumulation order, while
    # bounding harmonic/coordinate working memory for N=192..256 states.
    full_a,full_pa=np.zeros(n**3),np.zeros(n**3)
    for start in range(0,n**3,chunk_size):
        flat=np.arange(start,min(start+chunk_size,n**3))
        X,Y,Z=x[flat//(n*n)],x[(flat//n)%n],x[flat%n]
        r=np.sqrt(X*X+Y*Y+Z*Z);keep=r<(np.inf if r_cut is None else r_cut)
        if not keep.any():
            continue
        rk=r[keep]
        theta=np.arccos(np.clip(np.divide(Z[keep],rk,out=np.ones_like(rk),where=rk>0),-1,1))
        azimuth=np.arctan2(Y[keep],X[keep]);a,pa=np.zeros(len(rk)),np.zeros(len(rk))
        for l,entries in groups:
            harmonics=[real_sph_harm(l,m,theta,azimuth) for m in range(-l,l+1)]
            for mode,coefficients in entries:
                radial=np.interp(rk,mode['r'],mode['u']/mode['r'],right=0.)
                for y,(cq,cp) in zip(harmonics,coefficients):
                    a+=cq*radial*y;pa+=cp*radial*y
        full_a[flat[keep]]=a;full_pa[flat[keep]]=pa
    return [full_a.reshape(n,n,n),full_pa.reshape(n,n,n)]


def eigen_state(n, box, r0, wall_width, modes=None, amplitude=0., seed=0, r_cut=None, profile=None):
    s = np.zeros((n, n, n, 4))
    s[..., 0] = dilaton_sphere(n, box, r0, wall_width) if profile is None else profile(radius(n,box))
    if modes is not None and amplitude:
        a, pa = eigen_gas(modes, n, box, seed, r_cut)
        s[..., 1], s[..., 3] = amplitude*a, amplitude*pa
    return s


def laplacian(field, dx):
    return sum(np.roll(field, 1, ax)+np.roll(field, -1, ax)-2*field for ax in range(3))/dx**2


def reference(initial, cfg, times):
    """Independent semi-discrete DOP853 evolution on the same periodic stencil."""
    shape = initial.shape;dx = cfg['lSide']/shape[0]
    lam, dv, m2 = cfg['lambda'], cfg['delta_v'], cfg['mass_tv']**2

    def rhs(_, flat):
        s = flat.reshape(shape)
        fphi, fa = force(s[..., 0], s[..., 1], lam, dv, m2, cfg.get('mass_power', 2),cfg.get('vacuum_stiffening',0.))
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
    state = np.asarray(state, dtype='<f8')
    if state.ndim!=4:
        raise ValueError('Require a four-dimensional state')
    n = state.shape[0]
    if not 4<=n<=MAX_N or state.shape != (n, n, n, 4) or not np.isfinite(state).all() or not np.isfinite(box) or box <= 0:
        raise ValueError('Require finite cubic (N,N,N,4) state and positive box')
    with Path(path).open('wb',buffering=0) as f:
        f.write(b'PKIN0001'+struct.pack('<Id', n, box))
        state.tofile(f)


def _read(path, magic, cells):
    path=Path(path)
    with path.open('rb',buffering=0) as f:
        header=f.read(28)
        if len(header)!=28 or header[:8]!=magic:
            raise ValueError('Invalid magic/header')
        n,box,u=struct.unpack('<I2d',header[8:])
        if not 4<=n<=MAX_N or not np.isfinite([box,u]).all() or box<=0:
            raise ValueError('Invalid metadata')
        count=4*cells(n)
        if path.stat().st_size!=28+8*count:
            raise ValueError('Invalid byte count')
        data=np.fromfile(f,dtype='<f8',count=count)
        if len(data)!=count or f.read(1):
            raise ValueError('Invalid byte count')
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


def config(n, box, dt, duration, *, lam, dv, mass_tv, every=10, frames=0, slices=0, far_margin=6., mass_power=2, monitor_radius=0., vacuum_stiffening=0., tracking_inner=0., tracking_outer=0.):
    _check_power(mass_power)
    if not isinstance(n,(int,np.integer)) or not 4<=n<=MAX_N:
        raise ValueError(f'Require integer N=4..{MAX_N}')
    structured_wall.check_shape(vacuum_stiffening)
    if not np.isfinite(monitor_radius) or not 0<=monitor_radius<box/2:
        raise ValueError('Require 0<=monitor_radius<L/2')
    if (tracking_inner or tracking_outer) and not 0<tracking_inner<tracking_outer<box/2:
        raise ValueError('Require tracking radii both zero or 0<inner<outer<L/2')
    return dict(N=n, lSide=box, dt=dt, tMax=duration, expansion=False, evolver='PV2',
                **{'lambda': lam}, delta_v=dv, mass_tv=mass_tv, far_margin=far_margin,
                output_every=every, snapshot_every=frames, slice_every=slices,
                **({'mass_power': mass_power} if mass_power != 2 else {}),
                **({'monitor_radius': monitor_radius} if monitor_radius else {}),
                **({'vacuum_stiffening': vacuum_stiffening} if vacuum_stiffening else {}),
                **({'tracking_inner':tracking_inner,'tracking_outer':tracking_outer} if tracking_outer else {}))


def prepare_run(initial,cfg,root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    cfg = dict(cfg)
    write_initial(root/'input.bin', initial, cfg['lSide'])
    cfg.update(initial_state=str((root/'input.bin').resolve()), result_dir=str((root/'native').resolve()))
    (root/'input.json').write_text(json.dumps(cfg, indent=2)+'\n')
    text = '\n'.join(f'{k} = {str(v).lower() if isinstance(v, bool) else v}' for k, v in cfg.items())+'\n'
    (root/'input.in').write_text(text)
    return cfg


def run_prepared(root,executable=EXECUTABLE):
    root=Path(root);cfg=json.loads((root/'input.json').read_text())
    result = subprocess.run([str(executable), 'input='+str((root/'input.in').resolve())],
                            capture_output=True, text=True)
    (root/'run.log').write_text(result.stdout+result.stderr)
    if result.returncode:
        raise RuntimeError(f'Pocket solver failed: {root}/run.log\n{result.stderr[-1500:]}')
    return cfg


def run(initial,cfg,root,executable=EXECUTABLE):
    prepare_run(initial,cfg,root)
    return run_prepared(root,executable)


def initial_transfer_exact(input_path,snapshot_path,chunk_bytes=1024*1024):
    """Exact bounded-memory PKIN -> initial PKST byte comparison."""
    if not isinstance(chunk_bytes,int) or chunk_bytes<=0:
        raise ValueError('Require positive integer chunk_bytes')
    with Path(input_path).open('rb') as source,Path(snapshot_path).open('rb') as target:
        a,b=source.read(20),target.read(28)
        if len(a)!=20 or len(b)!=28 or a[:8]!=b'PKIN0001' or b[:8]!=b'PKST0001':
            return False
        n,box=struct.unpack('<Id',a[8:]);m,length,u=struct.unpack('<I2d',b[8:])
        if n!=m or box!=length or u!=0. or not 4<=n<=MAX_N:
            return False
        if Path(input_path).stat().st_size!=20+32*n**3 or Path(snapshot_path).stat().st_size!=28+32*n**3:
            return False
        while True:
            left=source.read(chunk_bytes);right=target.read(chunk_bytes)
            if left!=right:
                return False
            if not left:
                return True


def history(directory):
    return np.genfromtxt(Path(directory)/'native/history.csv', delimiter=',', names=True)
