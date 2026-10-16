# Stochastic nucleation checks

This workspace investigates stochastic bubble nucleation, false-vacuum
exclusion, periodic growth and arrival maps before coupling catalogues to
field evolution. The homogeneous and centered-bubble checks remain separate
reference cases. The first geometry-only ensembles are implemented in
[run_survival.py](run_survival.py), with analytic F(t) as the primary observable.
See [EXPERIMENT.md](EXPERIMENT.md) for the predeclared parameters and acceptance
criteria, and [results/README.md](results/README.md) for recorded findings.
The first stochastic catalogue driving CosmoLattice field evolution is now
implemented in [run_field.py](run_field.py), with parameters and validation
criteria in [FIELD_EXPERIMENT.md](FIELD_EXPERIMENT.md) and recorded
[field results](results/field/README.md). Bubble growth and expansion remain
prescribed, without axion backreaction.

Representative [bubble and axion GIFs](results/animations/README.md) use the
same frozen catalogue and an extended native field run. The reproducible
renderer is [make_gifs.py](make_gifs.py).
The [H*=0.05 GIFs](results/animations_H0.05/README.md) keep mt=1,
beta/H*=1 and H*L=3, with a newly evolved field; their detailed spatial
patterns are representative and have not passed a continuum-convergence check.
The [energy evolution](results/energy_H0.05/README.md) separates physical
kinetic, gradient and potential contributions with the paper's dimensionless
normalization and an independent homogeneous no-FOPT comparison.
The [Fig.2 physical-parameter rerun](results/paper_figure2/README.md) uses
H*/m=0.05,beta/H*=4,v=0.15,theta_i=2,L*m=8*pi,N=128 and sharp walls,
with fixed-catalogue timestep/grid checks. See its numerical-method caveats.

The first optimized engulfment implementation is
[spatial_engulfment.py](spatial_engulfment.py). Its periodic cell index replaces
only the frozen paper method's all-pairs tracer geometry. See the
[paired benchmark](benchmark_engulfment.py),
[predeclared contract](OPTIMIZATION_EXPERIMENT.md) and
[exact-equivalence/timing results](results/engulfment_index/README.md).
The original all-pairs path remains available for validation and fallback.

The warm-JIT tracer timing baseline is in
[benchmark_tracers.py](benchmark_tracers.py), with predeclared parameters in
[SCALING_EXPERIMENT.md](SCALING_EXPERIMENT.md) and measured phase/memory results
in [results/tracer_scaling/README.md](results/tracer_scaling/README.md).
It varies tracer population and box volume without changing the paper kernel.

```sh
MPLCONFIGDIR=/tmp/afopt-mpl python3 checks/stochastic/benchmark_tracers.py \
  --output outputs/stochastic/tracer-scaling-new
MPLCONFIGDIR=/tmp/afopt-mpl python3 checks/stochastic/benchmark_tracers.py \
  --density-match-only --output outputs/stochastic/tracer-density-new
```

The Hubble-rate scan is in [run_hubble.py](run_hubble.py), with its predeclared
contract in [HUBBLE_EXPERIMENT.md](HUBBLE_EXPERIMENT.md) and recorded results
in [results/hubble/README.md](results/hubble/README.md). It compares fixed
beta/H* scale covariance and fixed-beta changes to the transition shape.
Run with `MPLCONFIGDIR=/tmp/afopt-mpl python3 checks/stochastic/run_hubble.py
--output outputs/stochastic/hubble-new`, using a fresh output directory.

## Existing implementation

[src/axion_fopt_next/core.py](../../src/axion_fopt_next/core.py) already provides
seeded Poisson candidate sampling for a piecewise-linear conformal-time rate,
explicit rejection of births in converted space, periodic constant-speed bubble
growth and first-arrival maps. It uses brute-force geometry, which will remain
the correctness reference for acceleration.

Existing tests cover candidate-count mean/variance, linear-rate birth-time
distributions, exclusion, reproducibility, periodic causality and a constant-rate
false-vacuum ensemble. [workspace/statistics.json](../../workspace/statistics.json)
records that limited static-background check. It does not establish convergence
for an expanding radiation background or a physical transition parameter scan.

## Conventions we must preserve

Lengths are comoving, times are conformal, and c=1. The physical nucleation
rate Gamma is per physical volume per cosmic time. Candidate intensity is

    q(tau) = Gamma(t(tau)) a(tau)^4,
    E[N_candidates] = L^3 integral q(tau) dtau.

Sample candidates with q, then reject births in previously converted space.
Do not multiply q by the false-vacuum fraction as well: that would apply
exclusion twice. Candidate counts are Poisson; accepted counts generally are not.

Use one common constant physical wall speed initially:

    R_comoving = v_w (tau - tau_birth),
    tau_arrival(x) = min_j [tau_j + d_periodic(x,x_j)/v_w].

Growth begins at each bubble's own birth. Overlap denotes earliest conversion,
not a simulation of wall collision or plasma dynamics. Sampling positions and
births must remain independent of the axion lattice resolution.

For radiation domination, the repository conventions are

    a = H_* tau, t = H_* tau^2/2, t_* = 1/(2H_*),
    Gamma(t) = H_*^4 exp[beta(t-t_*)].

The sampler is exact for the supplied linear interpolant of q, not for an
arbitrarily coarse sampling of this exponential law. Rate-grid and time-window
convergence therefore need their own checks. See [physics.md](../../docs/physics.md).

## Order of investigation

1. **Candidate process.** Check integrated intensity, Poisson candidate counts,
   uniform positions, birth-time distributions, zero-rate intervals and resource
   guards. Separate RNG reproducibility from statistical agreement across seeds.
2. **Exclusion and periodic geometry.** Check causal accepted births, wraparound,
   first arrival and monotonic converted volume. For common speed, compare the
   accepted catalogue's conversion with all candidate bubbles (phantom
   equivalence). `sample_candidates` now exposes all attempted events, and
   `exclude_converted_births` applies the original brute-force rejection to those
   same events. `sample_catalogue` retains its existing seeded behavior.
3. **Static ensemble benchmark.** Revisit constant q with independent realizations
   and controlled geometry/probe error. In the regime R_max<L/2, compare against

       P_FV(tau) = exp[-I(tau)],
       I(tau) = (4*pi/3) v_w^3 integral q(s)(tau-s)^3 ds.

   Use ensemble uncertainty; neighboring lattice cells are correlated. Beyond
   R_max<L/2, replace Euclidean sphere volume with periodic torus-ball volume
   before comparing with an analytic survival law.
4. **Expanding rate benchmark.** Compare radiation-era survival and accepted
   comoving number density with independently integrated references. Check
   n_accepted(tau)=integral q(s)P_FV(s) ds using the survival law appropriate to
   the same box and initial time. Vary rate-grid spacing and start/end times;
   account for any omitted early nucleation tail.
5. **Catalogue and arrival-map convergence.** Separate ensemble size, box size,
   probe/lattice resolution and rate-grid resolution. Freeze catalogues when
   varying axion Ngrid. Define percolation/completion thresholds explicitly;
   do not assign a conventional number without stating it.
6. **Performance.** Measure sampling/exclusion and grid projection separately.
   Compare any spatial index against the same candidate events and brute-force
   arrival results. A shared seed alone does not guarantee identical events
   after changing rate-grid partitioning or RNG consumption.
7. **Field coupling.** After geometry/statistics are checked, evolve frozen
   catalogues with the validated adapter. Study a few seeded realizations before
   large ensembles. Test field timestep, lattice and wall-profile convergence
   independently of nucleation statistics.

The first new experiment should be a small geometry-only benchmark with a
fully specified rate, box, speed, time window and seed ensemble. Record its
acceptance criteria before examining its results. The first experiment uses a
static control and a radiation-era beta/H*=1 diagnostic; it is not the full
paper Fig. 1 scan.

## Run the initial comparison

Requires the project, SciPy, Matplotlib and Numba. The paper source is read from
the existing local checkout by default and verified against its pinned Git blob;
use `--paper-source` to provide that exact file elsewhere. Only four bubble and
tracer declarations are loaded, with JIT disk caching disabled. The original
field module is never imported and its files are not modified. Outer timesteps
are imposed by the geometry-only driver, rather than its field CFL controller.

```sh
PYTHONPATH=src MPLCONFIGDIR=/tmp/afopt-mpl python3 checks/stochastic/run_survival.py \
  --output outputs/stochastic/survival-new
```

Use a fresh output directory. Each case runs 24 realizations, two Poisson rate
grids and three tracer configurations. Both the internal tracer fraction and
independent geometric F are retained. The source snapshot, catalogues, seed
information, rates, reference curves, ensemble errors and checks are saved.

## Evidence and outputs

Each experiment should retain its exact rate table and physical parameters,
seed list/RNG identity and library versions, candidate/accepted counts,
catalogues, arrival-map conventions/checksums, independent references, measured
errors, ensemble uncertainty, numerical gates and runtime/memory measurements.
Use fresh directories under `outputs/stochastic/`; small recorded summaries,
figures and explanations can later live under `checks/stochastic/results/`.

Keep distinct: candidate versus accepted counts; probe discretization versus
seed-to-seed variation; free-growth radii versus collision-limited objects;
geometry validation versus field evolution; and a numerical check versus a
published abundance or spectrum reproduction.
