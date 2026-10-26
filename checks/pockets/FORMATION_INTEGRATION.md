# Formation potential: lightweight implementation checkpoint

2026-10-18. User requests limited computation while the laptop runs other work.
No lattice runs or full backend builds were launched for this checkpoint.

`formation_potential.h` implements a standalone scalar C++ kernel for the
potential defined in `paper_reference.py` and `PAPER_FORMATION_AUDIT.md`:

V = U(phi) + fa^2 mFV^2 exp(kappa phi) [1-cos(a/fa)].

Both fields have canonical kinetic terms in the chosen excursion units.
The tilted-quartic U is our explicit model choice. This is not a unique
potential specified by arXiv:2407.01676. Positive derivatives returned by
the kernel enter the equations as negative accelerations.

The kernel supplies separate vacuum and axion energies, both first
derivatives, the full symmetric Hessian, and axion-sector energy exchange.
For cosmic Minkowski time, allocating all V_ax to the axion sector gives

rho_a = (a_dot^2 + |grad a|^2)/2 + V_ax,
S_a = -a_dot grad a,
Q_a = (partial V_ax / partial phi) phi_dot,
partial_t rho_a + div S_a = Q_a.

The dilaton receives -Q_a. Q_a is internal exchange, not escape or total
energy loss. In FRW the dilution/work terms must be derived separately;
this checkpoint does not implement expansion.

## Small checks performed

`python3 -m unittest discover -s tests -p 'test_pocket_formation_cpp.py'`
passes. A single tiny C++ translation unit, compiled without optimization,
checks nine outputs against independent NumPy formulas at 69 configurations:
vacua, near-TV cancellation, finite FV backreaction, axion periodicity, and
64 fixed-seed samples. It also checks exponential overflow rejection and
invalid fa rejection. No framework or lattice evolution is involved.

`python3 -m unittest discover -s tests -p 'test_pocket_*reference.py'`
passes eight independent-reference tests, including finite-difference forces
and Hessians. These are model implementation checks, not physical survival
evidence. Tests use temporary files and leave no generated solver output.

## Homogeneous dynamics reference

The subsequent lightweight checkpoint adds `formation_homogeneous_probe.cpp`
and `validate_formation_homogeneous.py`. This is a standalone scalar C++ RK4
integrator, not CosmoLattice, and is not evidence for a spatial pocket.
It evolves both fields and integrates Q_a as a fifth ODE variable. An
independent Python/SciPy DOP853 solution supplies the comparison trajectory.

Fresh evidence: `outputs/pockets/formation-homogeneous-01/report.json`, raw
trajectories, build log, source and executable hashes. Fixed initial data
are (phi,a,phi_dot,a_dot)=(0,2,0,0.15), canonical excursion units, Minkowski
time 0..4. There is no expansion or spatial gradient in this check.

Predeclared gates were finest scaled trajectory error, relative total energy
drift and relative axion exchange residual each below 2e-7, and trajectory
error reduction above 10 for each timestep halving. All passed. For
steps 80/160/320, trajectory errors were 4.557e-7/2.931e-8/1.857e-9. The
finest relative energy drift was 3.612e-11; the axion exchange residual was
3.770e-11. Error scales are explicit in the report. Finite false-vacuum
axion energy drives the initially resting dilaton toward negative phi;
the trajectory is not artificially pinned at phi=0. This matters for future
initial-state preparation and is not a pocket-survival result.

This provides a checked reference for integration step 4 below, whose actual
compiled framework evolution remains pending. Run the validator with a new
`--output` directory to preserve existing evidence.

## Backend source integration checkpoint

The CPU backend now accepts `potential_model=formation` (default remains
`polynomial`). Formation selection requires explicit `mass_fv`, `mass_kappa`,
and `axion_fa`, and verifies `mass_tv=mass_fv*exp(mass_kappa/2)`. Polynomial
stiffening and nondefault `mass_power` are rejected with this model. These
are direct backend parameters; the existing Python production-run wrapper
has not yet been extended to prepare formation-model initial data.

Potential/force evaluation, the full mixed Hessian frequency diagnostic, and
both fixed and moving C++ energy monitors use the formation kernel when
selected. Independent NumPy monitor references accept a `formation` parameter
and pass instantaneous sector-continuity tests on a 12^3 inhomogeneous
state, including nonzero exchange and moving-boundary advection. Existing
polynomial reference-monitor tests also pass. These tiny algebraic checks
do not evolve a lattice.

A single CPU framework translation unit passed `-O0 -fsyntax-only`, using
the existing mapped-build include configuration, without replacing any
archived executable. Final-source evidence, including source hashes:
`outputs/pockets/formation-integration-syntax-02`.
The earlier syntax-01 attempted low-priority scheduling wrapper reported `setpriority: Operation
not permitted`; the compiler itself exited successfully. No new framework
executable, production build, or lattice run was produced. GPU execution is
not validated. Framework runtime checks and archived-run regression comparisons
remain required. The monitor-only runtime checkpoint below is now complete.

## Compiled monitor runtime checkpoint

To test the real monitor headers without compiling the full framework, the
unchanged scalar polynomial formulas and the optional formation dispatch
were moved from `pocket_model.h` into `pocket_potential.h`. This header is
shared by the actual model and both monitors. It carries the framework's
inline annotation when present and can also compile independently on CPU.
No substitute monitor or mock potential is used in this check.

`formation_monitor_probe.cpp` reads a 12^3 four-component field state and
executes the actual C++ fixed and moving monitors. The validator compares
six fixed outputs and twenty moving outputs with independent NumPy
references. It also evaluates the monitors at state +/- 1e-6 times the
independent semi-discrete field RHS, checking total and both sector
continuity equations, including moving-boundary advection and nonzero Q_a.
This directional-derivative check does not integrate a lattice trajectory.
The moving monitor's final FV enclosure-radius output is outside the twenty
reference quantities and is explicitly not counted as verified here.

Fresh successful evidence is in
`outputs/pockets/formation-monitor-agreement-02/report.json`, with the initial
state, nine raw C++ monitor outputs, compiler log and source/executable
hashes. The earlier attempt `formation-monitor-agreement-01` completed its
probes but failed JSON serialization of a NumPy boolean; it is retained as
an unsuccessful report-generation attempt, not a passing validation.

Predeclared normalized gates were 1e-12 for reference agreement and 1e-7
for instantaneous continuity. All six model/monitor pairs pass: formation,
polynomial phi^2, and polynomial phi^4 with vacuum stiffening 12, each in
fixed and moving windows. The largest reference discrepancy was 2.04e-16
of initial monitor energy; the largest continuity residual was 5.43e-11.
These are bookkeeping checks in dimensionless units, not time/spatial
convergence or evidence of long-term survival.

The scalar-kernel unit test now exercises actual potential dispatch orders
0/1/2 as well as the independent formation formulas. The extracted header
also passes the full CPU framework syntax check, with final source hashes
in `outputs/pockets/formation-integration-syntax-03`. No archived binary was
replaced and no framework evolution was launched.

## Actual framework homogeneous validation and polynomial regression

A single changed driver translation unit was compiled at O0 and linked to
existing CPU Kokkos/FFTW libraries, serially, in about two seconds. The new
executable is `build/pockets-formation-light-01/pocket_evolve`; its build
commands, source hashes and executable hash are in that directory's
`build.json`. Archived executables were not replaced. All validation runs
below explicitly used one OpenMP thread. This is a lightweight test build,
not a production build or evidence about performance.

`validate_formation_backend.py` ran the actual CosmoLattice coupled fields
on a homogeneous 4^3 grid, box 8, initial (phi,a,phi_dot,a_dot)=(0,2,0,0.15),
Minkowski cosmic time 0..4, with dt=0.02,0.01,0.005. Raw inputs, 21 snapshots
per run, histories, independent DOP853 trajectory, predeclared contract,
logs and hash provenance are in
`outputs/pockets/formation-backend-homogeneous-01`.

All declared checks passed. Scaled full-field trajectory errors were
2.6471e-5, 6.6160e-6 and 1.6539e-6: approximately fourfold reduction at each
timestep halving, as expected for PV2. Spatial spread was exactly zero.
The finest relative energy drift was 6.5977e-7; independent snapshot energy
agreed within 3.37e-16 of initial total energy. The frequency bound including
the mixed Hessian agreed within 6.94e-18 in absolute dt*frequency units.
Initial state transfer was exact. Six invalid configurations were rejected
before creating native output: inconsistent true-vacuum mass, zero fa,
unknown model, polynomial power/stiffening with formation, and missing mFV.
Their logs were inspected and identify the corresponding parameter errors.

The existing regression checker was then run at N=8, t=0..0.2, dt=0.005,
with a periodic-cut pocket against the archived mapped backend. All initial,
final, frame and slice binaries and all history bytes matched exactly for
p=2/4 and vacuum stiffening 0/12. Evidence:
`outputs/pockets/formation-polynomial-regression-01`. This is an implementation
regression check on a coarse small grid, not a new physical pocket result.
The four candidate runs took roughly 0.12 seconds each.

Neither homogeneous convergence nor coarse regression demonstrates spatial
convergence, survival, equilibrium, formation or cosmological scaling.
The next numerical gate is a small inhomogeneous formation-model comparison
with independent semi-discrete dynamics and integrated monitor continuity.

## Inhomogeneous framework dynamics and integrated accounting

`validate_formation_inhomogeneous.py` now exercises actual compiled coupled
fields at N=8, box=8, dx=1, t=0..0.2, dt=0.01/0.005/0.0025. It checks the
executable and all backend source hashes against the frozen light-build
manifest before running. All compiled runs explicitly use one OpenMP thread;
no new build or production simulation was needed.

The initial state is a smooth tanh dilaton pocket, Gaussian axion field and
momentum, plus translating and radial dilaton momentum. It is deliberately
a coarse diagnostic state, not a resolved scientific pocket. The field RHS
is independently integrated with DOP853 using the same periodic nearest-
neighbour spatial stencil. This holds discretisation fixed and tests time
integration, not spatial convergence. Initial amplitudes, state scales,
reference solution, parameters and contract are saved with every check.

Evidence: `outputs/pockets/formation-backend-inhomogeneous-01/report.json`,
all raw backend inputs/histories/frames, independent reference fields and
build provenance. All predeclared gates pass. Scaled trajectory errors are
7.7756e-5, 1.9438e-5 and 4.8595e-6. The finest relative total-energy drift
is 9.8030e-8. Independent snapshot total energy and fixed/moving monitor
outputs agree within 4.07e-16 of initial total energy.

Integrated total/dilaton/axion continuity includes outward currents, signed
exchange, and moving-boundary advection. The worst residual decreases from
1.6270e-6 through 4.0669e-7 to 1.0167e-7. Independently coarsening output
cadence on the finest trajectory changes integrated flux/exchange by at
most 2.36e-8 of initial monitor energy; this cadence audit is distinct from
timestep refinement. Nonzero moving-centroid speed (~0.175) and axion
exchange (~0.379) prevent a vacuous zero-current/exchange test.

An integer translation (3,-2,1) places the state across the periodic cut.
Final fields match the translated original exactly. Global total energy
and the invariant moving-monitor outputs agree within 2.23e-16 absolute.
The fixed monitor is intentionally not included in translation invariance,
since its physical location does not move with the state. The earlier
monitor-only checkpoint still supplies its independent continuity check.

These results establish a controlled implementation baseline for choosing
an accessible physical hierarchy. They do not establish wall or penetration
resolution, isolation, survival, equilibrium, or mass-radius scaling.

## Required before scientific field evolution

1. Runtime model selection/rejections and the small archived polynomial
   regression now pass, as recorded above. Retain these controls when changing
   initial-state preparation or later implementation.
2. Potential, derivatives and full mixed Hessian pass homogeneous validation;
   independent fixed-stencil inhomogeneous dynamics also now pass. Preserve
   this baseline while preparing physical initial data.
3. Compiled monitor/reference agreement, integrated continuity, independent
   output-cadence auditing, and periodic translation controls now pass on the
   small states above. Recheck these on each later physical hierarchy.
4. Homogeneous and inhomogeneous framework/ODE timestep-refinement checks now
   pass. Physical experiments must still refine space, time and boundaries
   separately at fixed physical parameters.
5. Predeclare an accessible scattering hierarchy, including wall-frame
   frequency, finite mFV, penetration and wall resolution. Only then run
   contracting fields and matched controls when compute is available.

The scalar kernel is connected to `PocketModel` and passes small homogeneous
and inhomogeneous framework validation, including integrated accounting.
Choosing and resolving a physical hierarchy is the next gate.
Existing toy-model simulation results remain results of their original
potential. No new trapping, radius, mass, stability or cosmology claim follows
from these implementation tests.
