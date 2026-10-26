# Single-bubble input parameters

For the later expanding single bubble and stochastic seeding pilot, use
[COSMOLOGY_INPUTS.md](COSMOLOGY_INPUTS.md). This file records the earlier static
experiments; its a=1/H=0 entries do not describe the cosmological extension.

This lists the entire physical and numerical setup for the longer experiment.
The later periodic-completion experiments have separate inputs and contracts:
[L24](SMALL_BOX_EXPERIMENT.md), [L20](SMALL_BOX_FOLLOWUP.md), and
[L20 refinements](SMALL_BOX_REFINEMENT.md). They keep the potential and seed
coefficients below but deliberately change the periodic geometry and duration.
The exact per-run native files and a machine-readable parameters.json are saved
under outputs/stochastic_v2/single-bubble-long-01/. See LONG_EXPERIMENT.md for
prospective acceptance criteria. This is an illustrative effective potential;
the parameter values are not inferred from the paper or matched to GeV scales.

## Units and fields

Natural units, c=1 (and hbar=1 for dimensional interpretation). The scalar
excursion sets the mass unit v=1. Fields S and canonical A have units v,
length/time units v^-1, potential densities units v^4, integrated energy units v.
theta=A/f is dimensionless. All tabulated values below use these units.
Restoring v, m_a^2=M^2(S/v)^2; the code uses v=1 explicitly. Kinetic terms
are (partial S)^2/2+(partial A)^2/2. There is no adjustable field-space metric.

## Physical model inputs

    U(S) = lambda/4 S^2(1-S)^2 + DeltaV(1-4S^3+3S^4)
    V(S,A) = U(S) + f^2 M^2 S^2 [1-cos(A/f)]
    S_tt = laplacian(S) - V_S
    A_tt = laplacian(A) - V_A

| Input | Value | Meaning |
|---|---:|---|
| v | 1 | Bare true-vacuum location / chosen field unit |
| lambda | 2 | Bare barrier-shape coupling |
| DeltaV | .06 | U(0)-U(1), true/false vacuum-energy difference |
| M | .5 | Axion mass at S=1 |
| f | .2 | Axion decay constant; A=f theta |
| a(t) | 1 | Static Minkowski background |
| H | 0 | No cosmic expansion |
| eta | 0 | No phenomenological friction |
| Axion force on S | included | Full reciprocal coupling |
| Plasma / thermal potential / noise | absent | No bath evolution or temperature input |
| Stochastic birth rate | absent | Exactly one initial bubble |
| Prescribed wall velocity v_w | absent | Velocity is measured from S evolution |
| Gravity / gravitational waves | absent | No metric or tensor evolution |

The bare vacuum energy zero is fixed: U(1)=0, U(0)=.06. There is no omitted
S-dependent offset. The degenerate-vacuum reference tension and wall-width
formula are diagnostics, not further independent inputs. The bare barrier is
at S=.367647, U=.0783863. Bare scalar masses are sqrt(U''(0))=1 and
sqrt(U''(1))=sqrt(1.72). With initial theta=1 the full true minimum is .994609;
the seed below deliberately uses the same profile as the short experiment.

## Initial condition inputs

One true-phase bubble centred at (L/2,L/2,L/2), with periodic nearest-image
radius r. No fields are inserted after t=0.

    S(r,0) = [1-tanh((r-R0)/w)]/2
    A(x,0) = f theta0
    S_dot(x,0) = A_dot(x,0) = 0.

| Input | Main / refinement / larger box | Empty control |
|---|---:|---:|
| R0 | 8 | 8 |
| Tanh width parameter w | 2 | 2 |
| theta0 | 1 radian | 0 |
| Canonical A0=f theta0 | .2 | 0 |
| S_dot0 | 0 | 0 |
| theta_dot0 | 0 | 0 |
| Bubble initial wall speed | 0 | 0 |
| Number of imposed bubbles | 1 | 1 |
| RNG seed | none; deterministic | none |

w=2 is not the 10%-90% thickness: the initial continuum thickness is
2*w*atanh(.8)=4.394449. The seed has tanh tails, not an exactly zero exterior
or an exactly homogeneous interior. No clipping is applied. It is not a
finite-temperature bounce, nor an exact static spherical solution.

## Numerical inputs

| Run | N per axis | Box L | dx | dt | t start / end | theta0 |
|---|---:|---:|---:|---:|---:|---:|
| Main | 128 | 48 | .375 | .04 | 0 / 18 | 1 |
| Coarser spatial control | 96 | 48 | .5 | .04 | 0 / 18 | 1 |
| Half-timestep control | 128 | 48 | .375 | .02 | 0 / 18 | 1 |
| Empty control | 128 | 48 | .375 | .04 | 0 / 18 | 0 |
| Larger-box control | 160 | 60 | .375 | .04 | 0 / 18 | 1 |

Each run starts from the same continuum formula, sampled at its cell centres;
the long run is not a restart from an interpolated old final state.

| Setting | Value |
|---|---|
| Spatial dimensions | 3 |
| Boundaries | periodic on all three axes |
| Grid locations | cell centres (i+.5)*L/N |
| Integration | actual CosmoLattice Position Verlet, PV2 |
| Laplacian | second-order nearest-neighbor finite difference |
| Energy gradient | compatible forward differences |
| Precision | double / float64 |
| Number of steps | 450 (dt=.04), 900 (dt=.02) |
| Histories | every .08: 2 or 4 integration steps |
| Central 2D slices | every .4: 10 or 20 steps; also final |
| Full 3D snapshots | initial and final |
| Background flags | expansion=false, fixedBackground=true, a0=1, t0=0 |
| Framework metadata | withGWs=false, omegaEoS=0, PS_type=1, kIR=2*pi/L |
| Backend | CPU, one process; existing OpenMP-capable libraries with one thread |
| Thread settings | OMP_NUM_THREADS=1, OMP_DYNAMIC=FALSE, OMP_PROC_BIND=false, OPENBLAS_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1 |
| Build | C++20, -O2, no fast-math; exact commands/hashes in build manifest |
| Per-run timeout | 180 seconds |

omegaEoS, PS_type and kIR are required/default framework metadata here: there
is no expansion, no spectrum estimator and no hydrodynamic fluid. Paths to
initial_state and result_dir are written as absolute native input parameters.
The driver requires fresh output directories and checks dt/dx<=.25 and
dt*sqrt(12/dx^2+largest nonnegative local potential-Hessian eigenvalue)<=.5
at diagnostic times; these are stability guards, not physical parameters.

## Measurement definitions (not evolution inputs)

Radius is the interpolated outward S=.5 crossing along a central coordinate
axis, correcting for the cell-centred transverse displacement. Thickness is
the distance between unique outward .9 and .1 crossings; ambiguous roots are
reported. Velocity is a finite-time derivative of measured radius, not a
terminal-speed model. No instantaneous velocity is supplied to the solver.
Energy is integrated over the box, with the scalar/axion split in the action.
Exchange is Q_A=2 f^2 M^2 S S_dot [1-cos(A/f)].

All initial conditions, potential inputs, boundary choices and time/grid inputs
must be held fixed in the comparisons except for the parameter named by the
control. Larger box changes its centre consistently and keeps the same physical
radial seed. No abundance, continuum limit or terminal wall speed is assumed.
