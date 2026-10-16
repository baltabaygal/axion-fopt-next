"""Independent Python references, initial states, and periodic defect diagnostics."""
from pathlib import Path
import json
import struct
import subprocess
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'build/pq_strings/pq_evolve'


def potential(x, y, lam, chi):
    r = np.hypot(x, y)
    return lam / 4 * (r*r-1)**2 + chi * (r*r-x*r)


def force(x, y, lam, chi):
    """Independent polar-coordinate gradient, with continuous zero origin limit."""
    r = np.hypot(x, y)
    c = np.divide(x, r, out=np.zeros_like(r), where=r > 0)
    s = np.divide(y, r, out=np.zeros_like(r), where=r > 0)
    radial = lam*(r*r-1)*r + 2*chi*r*(1-c)
    tangential = chi*r*s
    return np.stack((radial*c-tangential*s, radial*s+tangential*c), axis=-1)


def laplacian(field, dx):
    return sum(np.roll(field, 1, ax)+np.roll(field, -1, ax)-2*field
               for ax in range(3)) / dx**2


def coefficients(cfg, u):
    tau = cfg['tau_origin']+u
    a = tau/cfg['tau_origin'] if cfg['expansion'] else 1.
    if cfg.get('mass_schedule', 'powerlaw') == 'smooth_cap':
        s = np.clip((tau-cfg['mass_onset'])/cfg['mass_ramp_duration'], 0., 1.)
        chi = cfg['mass_cap']**2*s**3*(10.+s*(-15.+6.*s))
    else:
        chi = cfg['chi_star']*(tau/cfg['tau_star'])**cfg['n_tau']
    return tau, a, chi


def radial_coupling(cfg, u):
    tau = cfg['tau_origin']+u
    end = cfg.get('core_growth_until', 0.)
    physical = cfg.get('lambda', cfg.get('lambda_'))
    return physical*(end/tau)**4 if end > tau else physical


def reference(initial, cfg, times):
    shape = initial.shape
    dx = cfg['lSide']/shape[0]
    def rhs(u, flat):
        s = flat.reshape(shape)
        _, a, chi = coefficients(cfg, u)
        f, p = s[..., :2], s[..., 2:]
        acceleration = a*a*laplacian(f, dx)-a**4*force(f[..., 0], f[..., 1], radial_coupling(cfg, u), chi)
        return np.concatenate((p/a**2, acceleration), axis=-1).ravel()
    times = np.asarray(times, dtype=float)
    if np.any(times < 0) or np.any(times > cfg['tMax']):
        raise ValueError('Reference times outside evolution window')
    boundary = cfg.get('core_growth_until', 0.)-cfg['tau_origin']
    intervals = [0., cfg['tMax']]
    if 0 < boundary < cfg['tMax']:
        intervals.insert(1, boundary)
    values = np.empty((len(times), initial.size))
    state = initial.ravel()
    for start, end in zip(intervals[:-1], intervals[1:]):
        result = solve_ivp(rhs, (start, end), state, dense_output=True,
                           method='DOP853', rtol=2e-10, atol=2e-12, max_step=.02)
        if not result.success:
            raise RuntimeError(result.message)
        selected = (times >= start) & (times <= end)
        if selected.any():
            values[selected] = result.sol(times[selected]).T
        state = result.y[:, -1]
    return values.reshape((len(times),)+shape)


def write_initial(path, state, box):
    state = np.asarray(state, dtype='<f8')
    n = state.shape[0]
    if state.shape != (n, n, n, 4) or not np.isfinite(state).all() or box <= 0:
        raise ValueError('Require finite cubic (N,N,N,4) state and positive box')
    with Path(path).open('wb') as f:
        f.write(b'PQIN0001'+struct.pack('<Id', n, box))
        f.write(state.tobytes())


def read_state(path):
    raw = Path(path).read_bytes()
    if raw[:8] != b'PQST0001' or len(raw) < 52:
        raise ValueError('Invalid state magic/header')
    n, box, u, tau, a, chi = struct.unpack('<I5d', raw[8:52])
    if not 4 <= n <= 128 or not np.isfinite([box,u,tau,a,chi]).all() or box <= 0 or a <= 0 or tau <= 0 or chi < 0:
        raise ValueError('Invalid state metadata')
    if len(raw) != 52+32*n**3:
        raise ValueError('Invalid state byte count')
    state = np.frombuffer(raw, dtype='<f8', offset=52).reshape(n,n,n,4).copy()
    if not np.isfinite(state).all():
        raise ValueError('Nonfinite state')
    return dict(n=n, box=box, u=u, tau=tau, a=a, chi=chi), state


def config(n, box, dt, duration, *, lam=.08, origin=20., star=27., exponent=8.,
           chi_star=None, expansion=True, every=10, frames=0):
    return dict(N=n, lSide=box, dt=dt, tMax=duration, lambda_=lam,
                tau_origin=origin, tau_star=star, n_tau=exponent,
                chi_star=(origin/star**2)**2 if chi_star is None else chi_star,
                expansion=expansion, fixedBackground=True, H0=1/origin,
                omegaEoS=1/3, evolver='PV2', output_every=every, snapshot_every=frames)


def run(initial, cfg, root, executable=EXECUTABLE):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    cfg = dict(cfg)
    if 'lambda_' in cfg:
        cfg['lambda'] = cfg.pop('lambda_')
    write_initial(root/'input.bin', initial, cfg['lSide'])
    cfg.update(initial_state=str((root/'input.bin').resolve()), result_dir=str((root/'native').resolve()))
    (root/'input.json').write_text(json.dumps(cfg, indent=2)+'\n')
    text = '\n'.join(f'{k} = {str(v).lower() if isinstance(v,bool) else v}' for k,v in cfg.items())+'\n'
    (root/'input.in').write_text(text)
    result = subprocess.run([str(executable), 'input='+str((root/'input.in').resolve())],
                            capture_output=True, text=True)
    (root/'run.log').write_text(result.stdout+result.stderr)
    if result.returncode:
        raise RuntimeError(f'PQ solver failed: {root}/run.log\n{result.stderr[-1500:]}')
    return cfg


def energies(state, box, a, lam, chi):
    f, p = state[..., :2], state[..., 2:]
    kin = np.mean(np.sum(p*p, axis=-1))/(2*a**6)
    grad = sum(np.mean(np.sum((np.roll(f,-1,ax)-f)**2, axis=-1)) for ax in range(3))/(2*(box/len(f))**2*a*a)
    radial = np.mean(lam/4*(np.sum(f*f,axis=-1)-1)**2)
    tilt = np.mean(chi*(np.sum(f*f,axis=-1)-f[...,0]*np.linalg.norm(f,axis=-1)))
    return dict(kinetic=kin, gradient=grad, radial=radial, tilt=tilt, total=kin+grad+radial+tilt)


def wrap(angle):
    return (angle+np.pi)%(2*np.pi)-np.pi


def defects(state, box, tau):
    """Signed oriented plaquette windings and theta=pi crossings on links.

    2/3 orientation calibration follows the paper for isotropic networks;
    raw axis-aligned string lengths require separate interpretation.
    """
    if np.any(np.sum(state[..., :2]**2, axis=-1) == 0):
        raise ValueError('Phase is undefined at an exact zero: offset the diagnostic grid or use a core-aware interpolation')
    theta = np.arctan2(state[...,1], state[...,0])
    links = [wrap(np.roll(theta,-1,ax)-theta) for ax in range(3)]
    windings = []
    for i,j in ((0,1),(1,2),(2,0)):
        loop = links[i]+np.roll(links[j],-1,i)-np.roll(links[i],-1,j)-links[j]
        windings.append(np.rint(loop/(2*np.pi)).astype(np.int8))
    plaquettes = sum(int(np.abs(w).sum()) for w in windings)
    wall_links = sum(int((np.abs(np.roll(theta,-1,ax)-theta)>np.pi).sum()) for ax in range(3))
    dx = box/len(state)
    length = 2/3*plaquettes*dx
    area = 2/3*wall_links*dx*dx
    return dict(plaquettes=plaquettes, wall_links=wall_links,
                signed_winding=[int(w.sum()) for w in windings],
                string_length_estimate=length, wall_area_estimate=area,
                zeta=length*tau*tau/(4*box**3), A_dw=area*tau/(2*box**3)), windings


def deterministic(n, box, core=1., bias=0.):
    """Four neutral straight strings, periodic in all directions."""
    x = (np.arange(n)+.5)*box/n
    X,Y,Z = np.meshgrid(x,x,x,indexing='ij')
    # Z-independent analytic profile: two zeros of each factor -> four strings.
    real = np.cos(2*np.pi*X/box)+bias
    imag = np.sin(2*np.pi*Y/box)
    r = np.hypot(real,imag)
    scale = box/(2*np.pi*core)
    factor = np.divide(np.tanh(scale*r), r, out=np.full_like(r,scale), where=r>0)
    return np.stack((real*factor,imag*factor,np.zeros_like(r),np.zeros_like(r)),axis=-1)


def gaussian(n, box, correlation, seed):
    """Filter white noise by exp(-k^2*ell^2/4), giving continuum C(r)~exp(-r^2/(2ell^2))."""
    rng = np.random.default_rng(seed)
    k = 2*np.pi*np.fft.fftfreq(n, d=box/n)
    k2 = sum(a*a for a in np.meshgrid(k,k,k,indexing='ij'))
    filt = np.exp(-correlation**2*k2/4)
    fields = [np.fft.ifftn(np.fft.fftn(rng.normal(size=(n,n,n)))*filt).real for _ in range(2)]
    # Per-realization E[|Phi|^2]=1; each component has expected variance 1/2.
    f = np.stack(fields,axis=-1)
    f /= np.sqrt(np.mean(np.sum(f*f,axis=-1)))
    return np.concatenate((f,np.zeros_like(f)),axis=-1)


def cool(state, box, lam, duration=2., dt=.02):
    """Static chi=0 gradient flow; implicit lattice Laplacian, explicit potential.

    This pilot preparation is deliberately simpler than the paper's core-growth
    schedule. Cooling time is a fictitious coordinate, not conformal production time.
    """
    n = len(state);dx = box/n
    k = 2*np.pi*np.fft.fftfreq(n,d=dx)
    lattice_k2 = sum((2*np.sin(a*dx/2)/dx)**2 for a in np.meshgrid(k,k,k,indexing='ij'))
    steps = int(np.ceil(duration/dt));dt = duration/steps
    f = state[...,:2].copy();history = []
    for step in range(steps+1):
        s = np.concatenate((f,np.zeros_like(f)),axis=-1)
        e = energies(s,box,1.,lam,0.)['total']
        history.append((step*dt,e))
        if step<steps:
            nonlinear = lam*(np.sum(f*f,axis=-1)-1)[...,None]*f
            f = np.fft.ifftn(np.fft.fftn(f-dt*nonlinear,axes=(0,1,2))/(1+dt*lattice_k2[...,None]),axes=(0,1,2)).real
    return s,np.array(history)
