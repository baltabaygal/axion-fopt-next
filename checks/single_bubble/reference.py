"""Independent finite-difference reference and Fourier observables (m0=1)."""
from pathlib import Path
import struct
import numpy as np


def read_state(path):
    path = Path(path)
    with path.open("rb") as stream:
        if stream.read(8) != b"AFST0001":
            raise ValueError("Unknown state format")
        n, box, time, scale = struct.unpack("<Iddd", stream.read(28))
        if n < 1 or path.stat().st_size != 36 + 32*n**3:
            raise ValueError("Invalid state size")
        state = np.fromfile(stream, dtype="<f8").reshape(n,n,n,4)
    if not np.all(np.isfinite(state)):
        raise ValueError("Nonfinite single-bubble state")
    return dict(n=n, box=box, time=time, scale=scale), state


def mass(time, arrival, speed, width=.15):
    return .5 * (1 + np.tanh((time-arrival)*speed/width))


def integrate(arrival, box, speed, end, *, amplitude=1., dt_output=.02):
    """DOP853 of the periodic nearest-neighbor semi-discrete harmonic PDE."""
    from scipy.integrate import solve_ivp
    n = arrival.shape[0]
    count = arrival.size
    dx = box/n
    def rhs(time, state):
        field = state[:count].reshape(arrival.shape)
        laplacian = sum(np.roll(field,1,axis) + np.roll(field,-1,axis)-2*field
                        for axis in range(3))/dx**2
        return np.r_[state[count:], (laplacian-mass(time,arrival,speed)*field).ravel()]
    solution = solve_ivp(rhs,(0.,end),np.r_[np.full(count,amplitude),np.zeros(count)],
                         method="DOP853",rtol=2e-11,atol=2e-13,max_step=.01,dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    times = np.linspace(0,end,round(end/dt_output)+1)
    return times, solution.sol(times)


def observables(field, momentum, box):
    """Eq.55: sum (|pi_k|² + omega_k² |phi_k|²)/(2 omega_k), FFT/N³.

    Both continuum dispersion (paper) and stencil dispersion are returned.
    Complex modulus squared and ALL full-FFT modes, including k=0, are used.
    """
    n = field.shape[0]
    dx = box/n
    k = 2*np.pi*np.fft.fftfreq(n,d=dx)
    k2 = k[:,None,None]**2+k[None,:,None]**2+k[None,None,:]**2
    kd = (2*np.sin(k*dx/2)/dx)**2
    kd2 = kd[:,None,None]+kd[None,:,None]+kd[None,None,:]
    f2 = np.abs(np.fft.fftn(field)/field.size)**2
    p2 = np.abs(np.fft.fftn(momentum)/field.size)**2
    w = np.sqrt(1+k2)
    wd = np.sqrt(1+kd2)
    modes = (p2+w**2*f2)/(2*w)
    stencil_modes = (p2+wd**2*f2)/(2*wd)
    gradient = .5*sum(np.mean((np.roll(field,-1,a)-field)**2) for a in range(3))/dx**2
    energy_real = .5*np.mean(momentum**2+field**2)+gradient
    energy_fft = float(np.sum(.5*(p2+wd**2*f2)))
    edges = np.geomspace(2*np.pi/box*.8,np.sqrt(3)*np.pi/dx*1.001,40)
    shell, _ = np.histogram(np.sqrt(k2).ravel(),bins=edges,weights=modes.ravel())
    spectrum = np.column_stack([np.sqrt(edges[:-1]*edges[1:]),shell/np.log(edges[1:]/edges[:-1])])
    return dict(number_continuum=float(modes.sum()), number_stencil=float(stencil_modes.sum()),
                zero_mode_number=float(modes[0,0,0]), energy_real=float(energy_real), energy_fft=energy_fft,
                parseval_relative_error=float(abs(energy_real-energy_fft)/energy_real),
                mean_k=float(np.sum(modes*np.sqrt(k2))/modes.sum())), spectrum
