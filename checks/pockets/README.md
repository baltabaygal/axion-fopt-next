# Axion relic pockets

This workspace adapts the framework to look for axion relic pockets (ARPs):
false-vacuum regions in which axions are trapped because they are too light
to enter the surrounding true vacuum. The goal is to observe pocket-like
trapping in our lattice runs and compare with the reference paper where a
lattice can actually test it. A first dynamical-wall prototype is implemented
below.

The agreed field-only programme and the next predeclared pilot are in
[FIELD_GOAL_REVISED.md](FIELD_GOAL_REVISED.md), accepted on 2026-10-18.
It starts with the formation paper's model and hierarchy, then contracting
field pockets, survival/equilibration, scaling, expansion and collisions.
No particle or N-body implementation is in scope. The earlier static
programme remains in [FIELD_GOAL.md](FIELD_GOAL.md) as historical evidence.

Current field-only evidence (2026-10-19): finite-time stationary and moving
planar reflection benchmarks pass their declared controls; see
[PLANAR_STATIONARY_BENCHMARK.md](PLANAR_STATIONARY_BENCHMARK.md) and
[PLANAR_MOVING_BENCHMARK.md](PLANAR_MOVING_BENCHMARK.md). The moving result
uses the consistently validated optional fourth-order spatial operator;
earlier failures are retained. Neither is a demonstrated long-lived
three-dimensional pocket. Work remains bounded to one-thread planar tests
while the laptop is busy. A first finite-amplitude pilot observes sustained
wall reversal, now supported by temporal/spatial refinement and a displaced
second-wall control; see [PLANAR_PRESSURE_CONTROLLED.md](PLANAR_PRESSURE_CONTROLLED.md). Positive
vacuum bias now also retains the rebound under the same refinement and
boundary controls; see
[PLANAR_BIASED_PRESSURE_CONTROLLED.md](PLANAR_BIASED_PRESSURE_CONTROLLED.md).

Reference: W. Khater, M.C.D. Marsh and C. Nikolis, *The phenomenology of Axion
Relic Pockets*, [arXiv:2609.24714v1](https://arxiv.org/abs/2609.24714). The PDF
is retained here as [2609.24714.pdf](2609.24714.pdf) for lookups.
Formation theory is its Ref. [5]: P. Carenza, J. Eby, O. Iarygina and
M.C.D. Marsh, *Axion relic pockets — a theory of dark matter*, JHEP 09 (2024)
023, [arXiv:2407.01676](https://arxiv.org/abs/2407.01676) (not downloaded).

## What the paper contains

The paper is phenomenology. It has no lattice or field simulations; its
figures are a formation sketch, atomic fields, event-rate geometry and
detector limits.

Its formation reference, arXiv:2407.01676, does include a spherical
reflective-particle simulation (Sec. 5, Fig. 4). That is distinct from a
coupled-field simulation. We use its analytic relations for comparison;
the present programme uses fields throughout.

- **Sec. 2, Eqs. (2.1)-(2.7), from Ref. [5].** A dilaton first-order transition
  traps a pre-existing axion population in false vacuum. Bubble walls compress
  and heat the axions, which cannot enter the true vacuum. Their pressure
  stops the walls at the equilibrium radius R_eq. A pocket's energy density is
  4*DeltaV: 3*DeltaV in the hot axion gas and DeltaV of vacuum energy.
  Two parameters are left to simulations. For our cosmological convention use
  epsilon=a(t_transition)/a(t_matterlike), following Ref. [5], Eq. (16), with
  epsilon^4 <~ 0.15. The newer paper's prose instead references coalescence;
  these conventions must not be interchanged silently.
  alpha sets the wall Lorentz factor at coalescence, gamma_w = alpha*R_0/R_c,
  with 0 < alpha <~ 1.6.
- **Sec. 3, Eqs. (3.1)-(3.4).** Free axion in a totally reflective sphere:
  (d_t^2 - laplacian + m_a^2) a = 0 with a(R_pocket) = 0. Modes are
  j_l(k_nl r) Y_lm with k_nl = alpha_nl/R_pocket (alpha_nl the n-th zero of
  j_l) and omega^2 = k_nl^2 + m_a^2.
- **Sec. 3, Eqs. (3.28)-(3.36) and (3.47)-(3.48).** Classical and quantum
  axion-to-photon conversion for a pocket in a constant electric field E_0.
  The mode overlap integral gives F_nl(R) (3.30), single-mode power (3.31),
  summed power (3.36) and rate (3.47). In the k_a R >> 1 limit the rate
  becomes g^2 E_0^2 R/3 (3.48).
- **Secs. 4-6.** Atomic fields, IceCube, ANITA, Auger and direct-detection
  rates. These are outside lattice scope.

## Mapping to the existing code

The FOPT adapter already has the pocket mass structure. The axion is light
in false vacuum (`mass_false`) and heavy in true vacuum (`mass_true`). Arrival
maps accept `+inf` (never converted) and negative finite times (converted
before the start), so one static spherical false-vacuum pocket is a valid
input. Quadratic potential and `expansion=false` are supported.

Modes with omega < m_true are evanescent outside the pocket, with penetration
depth 1/sqrt(m_true^2 - omega^2). The paper's Dirichlet condition (3.2) is the
m_true -> infinity limit. Finite m_true, wall width and lattice spacing shift
the eigenfrequencies, and those shifts must be quantified.

Prescribed walls do not respond to axion pressure. Given a full arrival map,
false vacuum always disappears, so the existing model can show compression and
trapping but never stabilisation. Stabilised pockets need a dynamical wall:
a dilaton field whose motion is coupled to the axion mass.

## Prototype: dilaton+axion model with backreaction

A separate two-real-field CosmoLattice model makes the wall dynamical:

```text
V = lambda/4 phi^2 (1-phi)^2 + dV (1 - 4 phi^3 + 3 phi^4) + 1/2 m_TV^2 phi^2 a^2
```

Here phi is the dilaton (false vacuum 0, true vacuum 1) and a is a free axion,
massless in false vacuum. The axion pressure enters through the force
m_TV^2 phi a^2 on the dilaton. The first experiment evolves one ready-made
false-vacuum ball filled with a band-limited axion gas, on a static N=128
lattice, against no-gas and leaky (m_TV=0.3) controls.
The contract is in [EXPERIMENT.md](EXPERIMENT.md) and the results in
[results/README.md](results/README.md).

The default potential remains the above phi^2 coupling. An optional
`mass_power=4` replaces phi^2 by phi^4 in the axion mass, with independently
validated forces and energy accounting. This is an explicit alternative toy
model, not a claim to implement the paper's instanton mass function.
The [coupling contract](COUPLING_EXPERIMENT.md) and
[results](results/coupling-pilot-01.md) describe the comparison.

An optional `monitor_radius` (zero by default) adds a fixed-sphere energy
current and coupling-exchange diagnostic. It adds no damping or boundary
condition. Its symmetric-bond energy allocation is distinct from the old
forward-bond region energies; see [MONITOR_CONTRACT.md](MONITOR_CONTRACT.md).
Parameterized field suites can be run with `run_field_suite.py`; its
`--analyze-source` mode reconstructs reports in a fresh directory from
existing states and histories, without rerunning field evolution.
Optional `tracking_inner`/`tracking_outer` radii enable the validated moving
weighted-volume monitor; reports audit its corrected currents, cadence and
false-vacuum enclosure separately from the fixed monitor. See
[TRACKING_MONITOR_CONTRACT.md](TRACKING_MONITOR_CONTRACT.md) and the next
[global-stiffness experiment](GLOBAL_STIFFNESS_EXPERIMENT.md).
`initial_radial_balance.py` checks the seeded mode ensemble's finite-mass
radial force against independent radius differences, and gives a restricted
fixed-action, fixed-width reference radius. This is a mean spectral
diagnostic; field evolution must test stability and mass-radius scaling.
`compare_field_grids.py` also compares moving-window histories under the
[additional grid gates](MOVING_GRID_COMPARISON.md), when both histories
contain them. Prefix comparisons use shared endpoints and do not pass a
full-duration gate.

Outcome of prototype-01, in brief:
- Trapped axions hold the false vacuum open for more than four times the
  no-gas collapse time.
- The pocket is not stationary. About 20% of the axion energy escapes as
  axions upscattered above m_TV. The deformed wall radiates dilaton waves, and
  the pocket shrinks to R~11 by t=100. By the predeclared criteria, no relic
  pocket was established.

## Files

- `pocket_model.h`, `evolve.cpp`, `CMakeLists.txt`: compiled two-field model and PV2 driver, with strict guards.
- `model.py`: independent potential and force, initial states, gas modes, I/O, DOP853 and radial references, thin-wall ODE.
- `validate.py`: predeclared correctness gates for the compiled model.
- `run_prototype.py`: predeclared prototype suites, gates and outcomes.
- `plot.py`: figures and animation from saved runs.

## Build and run

From the repository root, using the cached pinned dependencies:

```sh
cmake -S checks/pockets -B build/pockets \
  -DCOSMOLATTICE_SOURCE="$PWD/build/vendor/cosmolattice" \
  -DFETCHCONTENT_SOURCE_DIR_TEMPLAT="$PWD/build/vendor/templat" \
  -DFETCHCONTENT_SOURCE_DIR_KOKKOS="$PWD/build/vendor/kokkos"
cmake --build build/pockets --target pocket_evolve -j 4
python3 checks/pockets/validate.py --output outputs/pockets/validation-new
python3 checks/pockets/run_prototype.py --output outputs/pockets/prototype-new
MPLCONFIGDIR=/tmp/pockets-mpl python3 checks/pockets/plot.py --source outputs/pockets/prototype-new \
  --control outputs/pockets/prototype-new --output outputs/pockets/plots-new
```

Output directories must be fresh. NumPy and SciPy are required; plotting also
needs Matplotlib and Pillow. Set `AFOPT_POCKET_EXECUTABLE` to include the
backend validation in `tests/test_pockets.py`.
