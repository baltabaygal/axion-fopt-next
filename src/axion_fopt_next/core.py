"""All positions are comoving; all times are conformal; c = 1.

Candidate intensity q(tau) = Gamma(t(tau))*a(tau)**4 per comoving volume.
The tabulated q is linearly interpolated: sampling is exact for that interpolant.
Only one common constant physical wall speed is supported in this first version.
"""
from dataclasses import dataclass
import numpy as np


def radiation_background(tau, h_star):
    """a=H_* tau, t=H_* tau**2/2, normalized to a(t_*)=1."""
    tau = np.asarray(tau, dtype=float)
    if not np.isfinite(h_star) or h_star <= 0 or np.any(~np.isfinite(tau)) or np.any(tau < 0):
        raise ValueError("Require H_*>0 and finite tau>=0")
    return h_star * tau, 0.5 * h_star * tau**2


@dataclass
class Catalogue:
    box_size: float
    wall_speed: float
    times: np.ndarray
    positions: np.ndarray
    candidate_count: int = 0

    def __post_init__(self):
        self.times = np.asarray(self.times, dtype=float).copy()
        self.positions = np.asarray(self.positions, dtype=float).copy()
        if not np.isfinite(self.box_size) or self.box_size <= 0:
            raise ValueError("box_size must be positive")
        if not np.isfinite(self.wall_speed) or not 0 < self.wall_speed <= 1:
            raise ValueError("Require 0 < wall_speed <= 1")
        if self.positions.shape != (len(self.times), 3) or self.times.ndim != 1:
            raise ValueError("Expected times (N,) and positions (N,3)")
        if np.any(~np.isfinite(self.times)) or np.any(~np.isfinite(self.positions)):
            raise ValueError("Catalogue must be finite")
        if np.any(np.diff(self.times) < 0):
            raise ValueError("Birth times must be ordered")
        self.positions %= self.box_size
        self.times.setflags(write=False)
        self.positions.setflags(write=False)

    def arrival_times(self, points, chunk_size=4096):
        """First wall arrival on a periodic cubic torus; empty catalogue => inf.

        Memory O(chunk_size), work O(N_bubbles*N_points). This is the reference
        geometry implementation, deliberately retained as an optimization oracle.
        """
        points = np.asarray(points, dtype=float)
        if points.ndim < 1 or points.shape[-1] != 3 or np.any(~np.isfinite(points)):
            raise ValueError("Expected finite points (...,3)")
        if not isinstance(chunk_size, int) or chunk_size < 1:
            raise ValueError("chunk_size must be a positive integer")
        flat = points.reshape(-1, 3) % self.box_size
        result = np.full(len(flat), np.inf)
        for start in range(0, len(flat), chunk_size):
            p = flat[start:start + chunk_size]
            arrivals = result[start:start + len(p)]
            for birth, centre in zip(self.times, self.positions):
                delta = np.abs(p - centre)
                delta = np.minimum(delta, self.box_size - delta)
                np.minimum(arrivals, birth + np.linalg.norm(delta, axis=1)/self.wall_speed, out=arrivals)
        return result.reshape(points.shape[:-1])

    def phase(self, tau, points):
        if not np.isfinite(tau):
            raise ValueError("tau must be finite")
        return self.arrival_times(points) <= tau

    def lattice_arrival_times(self, n):
        """Cell-centred grid, ndarray axes (x,y,z); no a**2 factor included."""
        if not isinstance(n, int) or n < 1:
            raise ValueError("n must be a positive integer")
        x = (np.arange(n) + 0.5) * self.box_size / n
        return self.arrival_times(np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1))


def sample_catalogue(tau, q, box_size, wall_speed, seed, max_candidates=1_000_000):
    """Generate Poisson candidates, then reject events in previously converted space.

    q is an attempted nucleation intensity, NOT Gamma multiplied by P_FV.
    Each interval has Poisson count V*integral(q); conditional birth times follow
    the normalized linear q. A resource guard prevents accidental huge runs.
    """
    tau, q = np.asarray(tau, float), np.asarray(q, float)
    Catalogue(box_size, wall_speed, np.empty(0), np.empty((0, 3)))
    if tau.ndim != 1 or len(tau) < 2 or q.shape != tau.shape:
        raise ValueError("tau and q must be matching 1D arrays with >=2 entries")
    if np.any(~np.isfinite(tau)) or np.any(np.diff(tau) <= 0) or np.any(~np.isfinite(q)) or np.any(q < 0):
        raise ValueError("Require increasing finite times and finite nonnegative q")
    dt = np.diff(tau)
    area = dt * (q[:-1]/2 + q[1:]/2)
    means = box_size**3 * area
    if not np.all(np.isfinite(means)) or means.sum() > max_candidates:
        raise ValueError("Expected candidates exceed budget; shorten interval or increase explicit budget")
    rng = np.random.default_rng(seed)
    counts = rng.poisson(means)
    if counts.sum() > max_candidates:
        raise ValueError("Realized candidates exceed budget")
    batches = []
    for i, count in enumerate(counts):
        if not count:
            continue
        y = rng.random(count) * area[i]
        slope = (q[i+1] - q[i]) / dt[i]
        # Rationalized quadratic root, valid for increasing or decreasing q.
        denominator = q[i] + np.sqrt(np.maximum(0, q[i]**2 + 2*slope*y))
        offset = np.divide(2*y, denominator, out=np.zeros_like(y), where=denominator > 0)
        batches.append(tau[i] + offset)
    births = np.sort(np.concatenate(batches)) if batches else np.empty(0)
    centres = rng.uniform(0, box_size, (len(births), 3))
    accepted_t, accepted_x = [], []
    for birth, centre in zip(births, centres):
        if accepted_t:
            delta = np.abs(np.asarray(accepted_x) - centre)
            delta = np.minimum(delta, box_size - delta)
            radii = wall_speed * (birth - np.asarray(accepted_t))
            if np.any(np.sum(delta**2, axis=1) <= radii**2):
                continue
        accepted_t.append(birth)
        accepted_x.append(centre)
    return Catalogue(box_size, wall_speed, np.asarray(accepted_t), np.asarray(accepted_x).reshape(-1, 3), len(births))


@dataclass(frozen=True)
class AxionPotential:
    """V/f**2=m**2(1-cos(theta)); force is d(V/f**2)/dtheta.

    The solver adds a**2 to this force in the conformal-time equation.
    This sharp mass step is a prescribed background, with external energy input.
    """
    mass_false: float = 0.0
    mass_true: float = 1.0

    def __post_init__(self):
        if not np.isfinite(self.mass_false) or not np.isfinite(self.mass_true) or min(self.mass_false, self.mass_true) < 0:
            raise ValueError("Masses must be finite and nonnegative")

    def mass_squared(self, tau, arrival):
        return np.where(tau >= np.asarray(arrival), self.mass_true**2, self.mass_false**2)

    def force(self, theta, tau, arrival):
        return self.mass_squared(tau, arrival) * np.sin(theta)

    def energy(self, theta, tau, arrival):
        return self.mass_squared(tau, arrival) * (1 - np.cos(theta))
