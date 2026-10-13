# CosmoLattice field-evolution adapter

Implemented and run: a single axion scalar evolved by CosmoLattice v2.0.0's
position-Verlet (PV2) integrator and TempLat's periodic lattice Laplacian.
The stochastic backend remains independent: the solver consumes an `AFOPT001`
arrival map, not nucleation code. An independent DOP853 reference checks the
homogeneous limits and a one-bubble solution at fixed spatial resolution.

This first adapter supports one CPU process, a fixed cubic 3D grid, homogeneous
initial conditions, and either no expansion or prescribed power-law expansion.
It does not reproduce a published axion abundance yet. MPI/GPU, restarts,
self-consistent expansion, other integrators, vacuum fluctuations and dynamical
transition fields are not implemented here.

## Build from the repository root

Dependencies: C++20 compiler, CMake >=3.20, Python >=3.10, Git, FFTW development
headers/libraries. CMake downloads TempLat v1.0.0 and Kokkos 5.1.1. Upstream
code is kept under ignored build storage, not committed into this repository.

```sh
mkdir -p build/vendor
git clone --depth 1 --branch v2.0.0 https://github.com/cosmolattice/cosmolattice.git build/vendor/cosmolattice
cmake -S adapters/cosmolattice -B build/cosmolattice \
  -DCOSMOLATTICE_SOURCE="$PWD/build/vendor/cosmolattice"
cmake --build build/cosmolattice --target axion_fopt_evolve -j 4
```

On the development machine, the checkout and dependencies are already cached
under `build/vendor`, and the executable is `build/cosmolattice/axion_fopt_evolve`.
For an offline reconfigure, also pass:

```sh
cmake -S adapters/cosmolattice -B build/cosmolattice \
  -DCOSMOLATTICE_SOURCE="$PWD/build/vendor/cosmolattice" \
  -DFETCHCONTENT_SOURCE_DIR_TEMPLAT="$PWD/build/vendor/templat" \
  -DFETCHCONTENT_SOURCE_DIR_KOKKOS="$PWD/build/vendor/kokkos"
```

The CMake console mentions lphi4 because upstream configures its default target;
our explicitly selected target compiles the AxionFopt model and driver.

## Run the deterministic bubble

```sh
PYTHONPATH=src python examples/export_one_bubble.py
build/cosmolattice/axion_fopt_evolve input=adapters/cosmolattice/one_bubble.in
```

The solver requires a fresh `result_dir` and refuses to overwrite existing
results. Override it on subsequent runs, for example
`result_dir=outputs/one_bubble_evolution_2`.

Results contain `history.csv` (time, scale factor, field and momentum averages,
mean potential), `initial.csv` and `final.csv` (every cell's field, canonical
momentum, arrival time and mass squared), and `parameters.txt`. CSV is a small
prototype output format, not a production snapshot format. For production
runs, replace it with scalable snapshots and observables.

## Units and initial data

Set alpha=1, fStar=f_axion. The evolved scalar is theta=phi/f_axion and
u=omega_star*(tau-tau_origin). The map uses conformal time and comoving length
in a common inverse-mass unit; `lSide=omega_star*map.box_size` and map n must
match N. The mass in the program potential is (m/omega_star)^2. The solver's
momentum is pi=a^2*dtheta/du. `theta_prime_initial` means dtheta/dtau; its initial
program momentum is theta_prime_initial/omega_star because a_initial=1.

Grid index (x,y,z) labels the physical point
(x+1/2,y+1/2,z+1/2)*L/N. This uniformly shifts the native periodic lattice origin
by half a cell without changing its stencil. Loading and snapshots use
TempLat's layout mapping, including its ghost padding. Validation compares every
transferred arrival cell and the final local mass to the original Python map.

For normalized radiation domination, choose tau_origin>0,
H0=1/tau_origin and omegaEoS=1/3. Then a=1+u/(omega_star*tau_origin)=tau/tau_origin.
The supplied bubble input uses tau_origin=omega_star=H0=1. That normalization
should not be confused with a physical H_* unless units have been matched.

## Timing hook and upstream pin

Official source: https://github.com/cosmolattice/cosmolattice/tree/v2.0.0
Pinned commit: `0dad16ec2c894a7c6abbea0e98734126b3470b36`.
TempLat v1.0.0: `f5239eebdfaa159e45b38b89d7ef8ef6e2222f4b`.

`prepare_overlay.py` verifies the upstream commit, clean tracked source and the
original position-Verlet header hash. It generates a header in the build tree,
leaving upstream files untouched, with three optional model hooks:

1. Set the external time immediately before a kick to tcurrent+w*dt/2.
2. Record tcurrent+w*dt after the final drift. The model defers the endpoint
   mass assignment until a potential/curvature is requested or measurements
   are synchronized; the next force midpoint can replace an unused endpoint.
3. Set it during measurement synchronization to the measured time.

The model updates a lattice mass field from the immutable arrival field when
needed. Potential, force and curvature all refresh a deferred mass before
returning an expression, and snapshots do the same. Identical time requests
reuse the mass field. The upstream
force and drift kernels supply the scale-factor powers. The driver restricts
execution to PV2; the generated hooks are not a validation claim for other
composition orders. Runge-Kutta and leapfrog require separate timing work.

TempLat adds `-ffast-math` in Release. This build removes it from TempLat's
interface options, since the arrival contract uses positive infinity and the
input guards depend on finite/NaN checks. The configured versions and backend
are recorded in `build-info.json` beside the executable.

The CPU adapter measures field moments and optional energies in one pass over
owned cells. It refreshes periodic ghosts before evaluating the same forward
gradient used by TempLat. The harmonic potential evaluates no unused sine or
cosine. Smooth-wall tanh calls are skipped only at |argument|>=20, where the
double-precision result is already rounded to +/-1; resolved wall values use
the original tanh formula. No float32 or timestep approximation is introduced.

For smooth profiles with at least eightfold repetition of arrival times, an
exact-value cache evaluates the mass once per distinct arrival value and
scatters it into the lattice. Equality uses the stored doubles without rounding
or radial binning. For the centered N=128 bubble, 14,167 distinct values replace
2,097,152 repeated tanh evaluations per stage. This adds roughly 9 MiB of index
and level storage. Set `mass_cache=false` to disable it; irregular maps with
too many distinct times automatically use the ordinary field expression.

Set `profile=true` to write `timings.json` with initialization, evolution,
diagnostic and snapshot durations. `mass_update_seconds` is a nested timing
included in the other phases and must not be added to their totals. For the
paired 128^3 benchmark, see checks/benchmark_runtime.py and
[recorded runtime results](../checks/runtime_results/README.md).

## Validation

SciPy is needed only for the independent reference validation:

```sh
python -m pip install -e '.[validation]'
PYTHONPATH=src python examples/validate_cosmolattice.py \
  --executable build/cosmolattice/axion_fopt_evolve
AFOPT_COSMOLATTICE_EXECUTABLE="$PWD/build/cosmolattice/axion_fopt_evolve" \
  PYTHONPATH=src python -m unittest discover -s tests -v
```

The validator requires a fresh output directory; pass `--output` to select one.
It checks zero-mass analytic drift, constant and delayed homogeneous masses,
a switch inside the first timestep, full arrival/mass transfers, scale factors,
endpoint potential measurements and a periodic one-bubble field. Reference
integration splits exactly at the mass discontinuities. Results are saved to
`validation.json`; [recorded development results](cosmolattice-validation.json)
include build identity and numerical errors.

Homogeneous constant and step-aligned turn-on cases show second-order temporal
convergence. Sharp unaligned bubble crossings generally reduce temporal order;
they must be resolved through smaller timesteps or a finite wall profile. The
one-bubble comparison holds N=8 fixed and is not a spatial continuum-convergence
claim. The driver bounds dt/dx but that alone does not guarantee resolution of
large masses; check temporal and spatial convergence for each physics setup.

Because the prescribed mass exchanges energy with the scalar, closed-system
energy conservation is not an acceptance criterion. Self-consistent Friedmann
feedback would require an explicit source-energy budget and is rejected.

## Plot the one-bubble run

For the non-expanding harmonic setup of Lee et al. (2024), with a smooth mass
wall and Fourier particle-number measurements, see
[checks/single_bubble](../checks/single_bubble/README.md). Optional adapter
parameters `potential=quadratic`, `wall_width>0` and `profile_wall_speed`
implement that equation. Defaults remain `potential=cosine` and `wall_width=0`
for existing runs. `energy_diagnostics=true` adds kinetic/gradient/total/inside
energy and external-power columns, and currently requires no expansion.

`snapshot_format=binary` selects AFST0001 state snapshots: 8-byte magic,
little-endian uint32 N, float64 program box length, program time and scale
factor (36-byte header), followed by N^3 interleaved float64 values
(theta,pi,arrival,mass2), C order with z fastest. Arrival values retain their
conformal units. `snapshot_at` optionally saves one checkpoint at the first
step at or after the requested program time. Existing CSV snapshots remain
the default. The binary writer currently supports little-endian hosts.

```sh
python -m pip install -e '.[plots]'
python examples/plot_cosmolattice.py --run outputs/one_bubble_evolution
```

Creates PNG and PDF versions of the mean field with spatial standard deviation,
and a bubble-plane slice of wall-arrival times and the final field. The shaded
band describes variation across cells, not uncertainty in the mean. The slice
uses the grid layer nearest the bubble's z coordinate. These are diagnostics of
the 16³ prototype run, not converged research figures.
