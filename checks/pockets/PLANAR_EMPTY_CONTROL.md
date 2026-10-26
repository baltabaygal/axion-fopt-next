# Short planar wall transient and diagnostic checks

2026-10-18. Actual compiled NDim=1 CosmoLattice coupled fields, one thread,
no axion packet. Evidence: `outputs/pockets/planar-empty-control-01/report.json`,
raw initial/final fields, histories, continuum reference final states, input
parameters, compiler/executable provenance and predeclared gates.

Four runs cover stationary and u=.3 inward-moving walls, at N=2048/4096,
L=512, dt=.025, t=0..4. Both fields evolve; exact zero axion data remain
zero. A strict 12-second total solver budget was declared; actual total
solver time was 2.30 seconds. The existing light executable was reused,
with its source and executable hashes checked before launching. These are
short transient controls, not full t=180/220 scattering controls.

All declared short-control gates pass. Stationary walls retain their
half-height position to numerical precision. Moving-wall position differs
from the continuum boosted kink-pair reference by .008179 on the coarse
grid and .002085 on the fine grid. Full-field errors fall from .005142 to
.001263; stationary full-field errors fall from .001646 to .000427.
This is evidence of improving spatial agreement at fixed timestep, not a
separate time-refinement result or proof of exact lattice stationarity.

Maximum relative energy drift over all four runs is 1.46e-7. Wall-width
changes stay below .23%. Axion field/momentum, axion energy and exchange
remain exactly zero in every control. Local phase-velocity estimates use
interpolated momentum divided by a bond gradient. Their interpolation bias
is visible even in the initial sampled continuum wall; they are contour
velocities, not independent material-velocity measurements. Position drift
and full-state comparison are therefore retained alongside that estimator.

These controls confirm that a transient exists and is modest over t<=4.
They do not bound the accumulated drift before packet arrival at t~100/131.
Matched empty-wall evolution over the actual measurement interval and
independent timestep refinement remain necessary before claiming weak-axion
backreaction or a quantitative scattering coefficient. A long-window cost
extrapolation from these short runs is only a planning estimate; no further
runs were launched for this checkpoint.

## Independent geometry and local energy accounting

`planar_diagnostics.py` locates the unique rising phi=.1,.5,.9 crossings,
unwraps them locally, and rejects ambiguous or incorrectly ordered widths.
It does not silently choose one of multiple rising walls. Tests compare a
boosted continuum seed's width and velocity and exercise ambiguity rejection.

Its independent `window_budget` splits nearest-neighbour bond gradient
energy equally between endpoints. For weights W_i and derivatives dW_i/dt,
the outgoing lab flux is sum_i J_(i,i+1)(W_i-W_(i+1)), with
J=-(f_(i+1)-f_i)(pi_i+pi_(i+1))/(2dx). Boundary advection is
sum_i (dW_i/dt) rho_i dx. Corrected flux is lab flux minus advection.
The weighted axion source is sum_i W_i kappa V_ax pi_phi dx.

Finite-difference differentiation along an independent semi-discrete field
RHS checks total, dilaton and axion moving-window continuity, with explicitly
nonzero exchange. Both geometry/accounting tests pass in milliseconds.
This provides a checked independent diagnostic. It has not yet been wired
into the compiled line driver's output, whose histories still contain only
global energies and exchange. Scattering fractions are not yet measured.

The next step is to expose matching local diagnostics during the actual
compiled evolution and run one bounded weak-packet comparison. This planar
work supports local reflection tests; it does not establish three-dimensional
pocket survival, scaling, expansion or formation.
