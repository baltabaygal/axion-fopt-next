# Modest-grid PQ strings and domain walls

This workspace develops a small post-inflationary complex-PQ simulation:
formation/evolution of global strings followed by N_dw=1 domain walls and
collapse as the axion mass grows. The first goal is validated, observable
defect dynamics on a laptop-sized periodic grid. A separate compiled two-field
CosmoLattice model, independent references, and 32³/64³ experiments are now
implemented. See [recorded results](results/README.md).

Reference: Correia et al., *The axion mass in the post-inflationary scenario
in a minimal scaling model*, [arXiv:2609.20801v1](https://arxiv.org/abs/2609.20801),
especially Eqs. (1)-(4), (45), (47)-(48), and Sec. IV.
The downloaded [paper](../../papers/2609.20801.pdf) is retained in this repo.
This experiment does not aim to reproduce its production-scale axion mass
prediction or its statistical precision.

## Physical model and conventions

Use Phi = X + iY, r = sqrt(X^2+Y^2), theta = arg(Phi), with kinetic term
1/2*((partial X)^2+(partial Y)^2). The paper's potential is

```text
V = lambda/4 * (r^2 - f_a^2)^2
    + chi(tau) * r^2/f_a^2 * (1 - cos(theta)).
```

The tilt in Cartesian coordinates is chi/f_a^2*(r^2-X*r), so its force can
be evaluated without a phase branch cut. Its gradient has a zero limit at
r=0; the implementation must handle that point explicitly. Near the
untilted vacuum the radial mass is m_s^2=2*lambda*f_a^2 and the angular
mass is m_a^2=chi/f_a^2. The tilt also changes radial dynamics away from
theta=0; fixing r=f_a would remove required physics.

For each component F in {X,Y}, use conformal time tau and comoving x:

```text
F'' + 2*(a'/a)*F' - laplacian(F) + a^2*dV/dF = 0.
pi_F = a^2*F'
F' = pi_F/a^2
pi_F' = a^2*laplacian(F) - a^4*dV/dF
```

Choose f_a=1 and prescribed radiation expansion a=tau/tau_ref with
tau_ref>0. Physical H=a'/a^2=tau_ref/tau^2. Set
chi(tau)=f_a^2*H(tau_star)^2*(tau/tau_star)^n_tau, which enforces
m_a(tau_star)=H(tau_star); n_tau=8 is used in the recorded network demonstration.
Physical lambda stays constant during production evolution. Preparation
stages may alter the equations only through an explicitly documented schedule.
No stochastic FOPT bubbles, arrival-map mass switches or gravitational
backreaction are used in this model.

## Backend and size

Reuse the pinned CosmoLattice v2.0.0 / TempLat infrastructure, periodic
nearest-neighbor Laplacian and force-stage timing approach. Add a separate
two-real-field model; the existing axion FOPT model remains independently usable.
Begin with CPU/PV2. The reference paper uses leapfrog, so matching its
integrator is not assumed.

The recorded pilot uses N=32, L=12; refinement uses N=64 at the same L.
The main network demonstration uses N=64, L=24, lambda=0.08, tau_origin=20,
tau_star=27, n_tau=8, tau_end=32 and dt=0.02 (with a dt=0.01 repeat).
These are deliberately reduced parameters chosen together using resolution
constraints, rather than the paper's production parameters.
Aim for at least four cells per characteristic comoving core width
1/(a*m_s) and wall width 1/(a*m_a) while defects are being measured; refine
to check whether this is sufficient. The box must also contain the intended
defects and avoid light-crossing contamination during the comparison window.
These constraints may limit the duration of a small-box run.

Simply retaining the paper's times, dx and correlation lengths while shrinking
its box would destroy its scale hierarchy. Our setup will document reduced
dynamic range and any modified preparation parameters explicitly.

## Implementation and validation sequence

1. Derive Cartesian potential, force and curvature; check force derivatives
   away from the origin, the origin limit, vacuum and small-amplitude masses.
2. Build the separate compiled model with homogeneous and arbitrary spatial
   initial-state loading. Validate homogeneous radial/angular evolution against
   an independent ODE, including time-dependent chi and expansion.
3. Compare a small periodic lattice against an independent implementation of
   the same semi-discrete equations. Check transfer, time mapping, refinement
   and energy/source accounting before evolving defects.
4. Prepare a deterministic periodic string–antistring pair or closed loop,
   first at chi=0, then with growing chi.
   Verify winding detection and inspect resolved core/wall evolution.
5. Add seeded correlated Gaussian complex initial fields and diffusion cooling.
   Investigate the paper's referenced core-growth preparation before adopting
   it. A direct cooled-network pilot may precede it, labelled as a simplification.
6. Record network evolution on N=32 and N=64 with independent temporal checks;
   only then add current-based spectra and exploratory particle-number estimates.

Before each run, save a concrete experiment contract with parameters,
initial-state definition, tolerances, observables and stopping criteria.
Neither defect disappearance alone nor one spectrum establishes cosmological
scaling or a relic abundance.

## Outputs and provenance

Raw runs will use fresh outputs/pq_strings/<run-name>/ directories, containing
inputs, initial states, snapshots, histories, logs, build/source hashes and
validation.json. Compact figures, tables and summaries can be retained under
checks/pq_strings/results/. Record phase preparation separately from physical
evolution and distinguish visual evidence from numerical acceptance gates.

The deterministic defect checks, annihilating configuration, and short random
network demonstration are complete. The observed loss of strings also occurs
in the massless control, so this short run does not isolate the paper's
wall-driven late collapse or establish network scaling.

## Build and run

From the repository root, using the already cached, pinned dependencies:

```sh
cmake -S checks/pq_strings -B build/pq_strings \
  -DCOSMOLATTICE_SOURCE="$PWD/build/vendor/cosmolattice" \
  -DFETCHCONTENT_SOURCE_DIR_TEMPLAT="$PWD/build/vendor/templat" \
  -DFETCHCONTENT_SOURCE_DIR_KOKKOS="$PWD/build/vendor/kokkos"
cmake --build build/pq_strings --target pq_evolve -j 4
python3 checks/pq_strings/validate.py --output outputs/pq_strings/validation-new
python3 checks/pq_strings/run_demo.py --suite pilot --output outputs/pq_strings/pilot-new
python3 checks/pq_strings/run_demo.py --suite refinement --output outputs/pq_strings/refinement-new
python3 checks/pq_strings/run_demo.py --suite matched --output outputs/pq_strings/matched-new
python3 checks/pq_strings/run_demo.py --suite collapse --output outputs/pq_strings/collapse-new
python3 checks/pq_strings/run_demo.py --suite network --output outputs/pq_strings/network-new
MPLCONFIGDIR=/tmp/pq-mpl python3 checks/pq_strings/plot.py \
  --source outputs/pq_strings/network-new --output outputs/pq_strings/plots-new
```

All output directories must be fresh. Each experiment writes its contract before
running and exits nonzero on a failed numerical gate. Requires NumPy and SciPy;
plotting additionally requires Matplotlib and Pillow. Existing CosmoLattice
compiler/FFTW requirements are described in [the adapter guide](../../docs/cosmolattice.md).
No new package dependency is imposed on the solver-independent core.

Enable both optional backend test suites:

```sh
AFOPT_COSMOLATTICE_EXECUTABLE="$PWD/build/cosmolattice/axion_fopt_evolve" \
AFOPT_PQ_EXECUTABLE="$PWD/build/pq_strings/pq_evolve" \
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

The `matched` suite prepares one N=64 random initial state and gives N=32
its 2×2×2 cell averages, holding L, physics and dt fixed. Separate random
seeds at different grids would not constitute matched initial data.
The `collapse` suite offsets the real part of the analytic four-string
configuration by +0.5, making opposite-string separations unequal. Both
this configuration and its massless control annihilate; the symmetric
configuration remains a useful symmetric topology/transfer check.

## Files and numerical contracts

- `pq_model.h`: two-field potential expressions and correctly staged chi(tau).
- `evolve.cpp`: compiled PV2 driver, strict input guards, diagnostics and snapshots.
- `model.py`: independent polar-force/ODE reference, initial-state preparation,
  periodic winding and branch-cut diagnostics, binary I/O.
- `validate.py`: independent numerical gates, including exact-zero force behavior,
  both signs of the Cartesian components, and expanding/inhomogeneous evolution.
- `run_demo.py`: frozen-parameter experiment suites, refinement, controls and provenance.
- `plot.py`: phase/modulus slices, diagnostics and an animation from saved solver states.

At r>0 the diagonal potential curvatures are

```text
V_XX = lambda*(3X^2+Y^2-1) + chi*(2-3X/r+X^3/r^3)
V_YY = lambda*(X^2+3Y^2-1) + chi*(2-X/r+X*Y^2/r^3)
V_XY = 2*lambda*X*Y - chi*Y^3/r^3.
```

The tilt is continuously differentiable at the origin but has no unique Hessian
there. PV2 only requires the finite gradient; no dummy origin curvature or
vacuum-fluctuation initialization is supplied. The driver checks a conservative
field-dependent frequency bound as well as dt/dx. Neither bound substitutes
for timestep refinement. Topology diagnostics reject exact-zero grid vertices,
where phase is undefined, rather than assigning an artificial winding.

Initial state `PQIN0001`: 8-byte magic, little-endian uint32 N, float64 box
length, then N³ interleaved float64 (X,Y,pi_X,pi_Y), C order with z fastest.
Header length is 20 bytes. State snapshot `PQST0001` adds float64
(u,tau,a,chi) after the box length, for a 52-byte header; component data
have the same order. Coordinates are cell-centered and all runs use
f_a=omega_star=1. Snapshots are synchronized at their stored time.
Loading a state does not implement a general restart with arbitrary a_initial.

Physical energy density is K+G+V, where K=(pi_X²+pi_Y²)/(2a⁶) and
G is the periodic forward-difference gradient energy divided by a². The
independent energy balance is

```text
rho' = -(a'/a)*(6K+2G) + chi'(tau)*<r^2-X*r>.
```

Static chi=0 diffusion uses implicit lattice Laplacian and explicit quartic
force, with fictitious dt=0.02 for time 2; its energy decreases monotonically.
It is a preparation simplification, not the reference paper's full cooling and
core-growth procedure. The Gaussian filter is exp(-k² ell²/4), ell=2.5;
before cooling each realization is normalized to <|Phi|²>=1, with zero momentum.

Raw wall-link counts locate the theta=pi branch surface, including in chi=0
controls where this surface is not a physical domain wall. The paper's 2/3
orientation factors are provided as exploratory isotropic-network estimators;
they underestimate a known axis-aligned straight string. Report raw counts
alongside these estimators. Phase/modulus movies are diagnostics, not an
axion-particle-number or abundance measurement.


## Yield audit and the capped-mass toy

The next phase implemented `spectra.py`: modulus-weighted angular currents,
midpoint forward links, aligned lattice-momentum longitudinal projection,
continuum and lattice dispersion estimates, and exact finite-mode shell sums.
The homogeneous mode is explicitly excluded from the current-only estimator.
Independent analytic modes and a compiled, evolving small-amplitude wave validate
normalization, Parseval identities and late adiabatic behavior. The short network
runs remain transient diagnostics; they do not supply a late abundance.

`feasibility.py` audits simultaneous core, wall, mass-hierarchy, horizon and IR
requirements. With the paper's growing mass law and tau/tau_star=3, four cells
per core and m_a/m_s<=0.16 require at least N=36,450 for the stated necessary
bounds. This is a fixed-grid necessary bound, not a hardware forecast or a proof
that another method cannot work. Core growth helps preparation but does not
remove the late physical scale hierarchy. `prepare_network.py` validates a
continuous lambda∝a^-4 preparation handoff and shows that the larger L=24
massless preparation retains a network; the L=12 example does not.

The user therefore selected an explicitly modified small-grid toy. It retains
the two-field potential and radiation expansion and replaces the power-law
susceptibility with a smooth transition to a constant physical mass:

```text
s = clip((tau - mass_onset)/mass_ramp_duration, 0, 1)
chi = mass_cap^2 * (10*s^3 - 15*s^4 + 6*s^5)
chi' = mass_cap^2 * 30*s^2*(s-1)^2 / mass_ramp_duration
```

The derivative vanishes outside the transition. The default `powerlaw` schedule
continues to reproduce the earlier experiments. `validate_cap.py` checks the
modified force stages and explicit energy source against an independent
semi-discrete ODE, including unaligned transition times.

The first toy uses L=24, N=64, lambda=0.08, tau_origin=40, core growth ending
at tau=44, a six-unit mass turn-on, and mass_cap=0.06. Thus m_a/m_s<=0.15.
Static diffusion prepares four symmetric neutral straight strings. A massless
control uses the exact same frozen state. N=128 uses Fourier interpolation of
that N=64 state with all physical parameters fixed; dt refinement is separate.
Seeded Gaussian networks are separate finite-volume realizations, not a scaling
ensemble. Raw runs and their predeclared gates are under outputs/pq_strings/toy-*.

`run_toy.py` records the contract and provenance before evolution. `analyze_toy.py`
reports topology, core/wall widths, energy components and current spectra.
For a defect-free field with max|theta|<=0.3 and radial RMS displacement<=0.1,
`coherent_zero_mode` separately adds the harmonic energy of the mean angle and
mean temporal current. It never infers the zero-mode potential from spatial
currents. Ineligible states have no completed yield.

The toy normalization is nu_toy=(a/a_star)^3*n_a/H_star, f_a=1, with tau_star
solving m_a=H during the capped turn-on. Equal dimensionless notation does not
make this toy a paper reproduction or a cosmological relic-density prediction.
A defensible toy yield requires two late mass periods, a <=5% plateau variation,
<=1% temporal and <=10% spatial change, and resolved defects through collapse.
Raw PQ total energy includes radial excitations and must be distinguished from
the harmonic angular energy proxy.

Example (fresh output directories required):

```sh
python3 checks/pq_strings/validate_cap.py --output outputs/pq_strings/cap-validation-new
python3 checks/pq_strings/run_toy.py --output outputs/pq_strings/toy-new --end 400
python3 checks/pq_strings/analyze_toy.py outputs/pq_strings/toy-new
```
