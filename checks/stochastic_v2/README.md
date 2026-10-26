# Dynamical transition field: research workspace

The 2026-10-20 request for longer stochastic evolution is recorded in
[LONG_STOCHASTIC.md](LONG_STOCHASTIC.md). The
[completed long-run assessment](results/long-stochastic-01/README.md) extends
to u=80,a=3 with grid and time controls. All cells on the main grid were
observed in the true phase at least once, but no simultaneous all-cell interval
exists in any run. Local recrossings remain and their small volume is grid
dependent; the persistent-conversion goal is unmet. The
[intermediate u40 checkpoint](results/long-stochastic-40-01/README.md) reaches
99.57% S>=.5, with no all-cell interval; this is not completed conversion.
Later stages and their controls must be assessed separately.

The requested expansion and stochastic-seeding extension is implemented and
validated. Start with the [cosmology inputs](COSMOLOGY_INPUTS.md),
[expanding single-bubble findings](results/expanding-single-bubble-01/README.md),
and [stochastic-pilot findings](results/expanding-stochastic-01/README.md).
The [prospective contract](EXPANSION_STOCHASTIC_CONTRACT.md) distinguishes
prescribed radiation expansion and finite-radius Poisson seeding from a
thermal/bounce calculation. The original polynomial coupling and reciprocal
axion force are retained. There is no imposed wall speed. Earlier research
notes and static evidence below remain historical context.

Prepared 2026-10-19 as a research workspace. The user subsequently authorized
the first single-bubble experiment. It is now implemented and its declared
pilot checks pass: see [the contract](EXPERIMENT.md) and
[recorded findings](results/single-bubble-02/README.md). This folder is a sibling
of `checks/stochastic`, not its replacement. The research discussion below is
retained; the experiment contract identifies the specific model actually run.

The requested longer evolution now reaches t=18. See the
[long-run findings](results/single-bubble-long-01/README.md),
[complete input parameters](INPUT_PARAMETERS.md), and
[long-run contract](LONG_EXPERIMENT.md). The original energy-summation failure
is preserved; compensated endpoint replay repairs that diagnostic without
replacing the long trajectories.

The [energy-evolution and legacy comparison plots](results/energy-comparison-02/README.md)
use the saved long-run histories and the existing `checks/single_bubble` runs.
They include explicit normalization and model differences; no matched-control
simulation or causal enhancement/suppression claim is implied by the overlay.

The later [smaller periodic-box experiments](results/periodic-conversion-01/README.md)
test filling the box. L20 reaches an all-cell S>=.5 interval for about one
time unit, confirmed by spatial/time refinements, then recrosses. This is
transient conversion, not a settled uniform vacuum; the original sustained
completion goal and failed pilot gates remain recorded.

## The question

Can we replace the prescribed bubble background with a dynamical scalar order
parameter S, and determine when that changes axion production?

There already is an evolving field in the old calculation: the axion
phi = f theta. The missing dynamical degree of freedom is S, which makes the
transition happen and controls the axion mass. Conceptually:

    current: sampled births -> prescribed expanding spheres -> m_a^2(x,t) -> theta
    proposed: seeded bubbles -> evolve S <-> evolve theta, with m_a^2 = m_a^2(S)

The arrows become two-way only if we retain the axion force on S. Evolving S
while omitting that force is a spectator approximation, not a fully coupled
model. It can be useful when its smallness has been established.

## Paper and source map

The pinned local paper is [2607.01333v2.pdf](references/2607.01333v2.pdf),
with [searchable text](references/2607.01333v2.txt) and
[download provenance](references/README.md).

Relevant paper locations: Sec. II for prescribed nucleation and growth;
Sec. III, Eqs. (4)-(10), for the S-dependent axion potential and its
thin-wall replacement; Appendix A, Eqs. (A1)-(A9), for an illustrative
renormalizable origin. Appendix A supplies motivation, not a complete
finite-temperature potential or transport calculation.

Read these existing repository sources before future implementation:

| Source | Role and limitation |
|---|---|
| [stochastic README](../stochastic/README.md) and [physics contract](../../docs/physics.md) | Existing conventions and scope |
| [core.py](../../src/axion_fopt_next/core.py) | Poisson candidate process, explicit exclusion, constant-speed geometry, sharp axion mass |
| [tracer_reference.py](../stochastic/tracer_reference.py) | Separate pinned paper tracer implementation; finite tracer and time-step effects |
| [run_field.py](../stochastic/run_field.py) | Frozen catalogue driving a cosine axion, plus independent semi-discrete reference |
| [axion_fopt.h](../../adapters/cosmolattice/axion_fopt.h) | Actual external-arrival-map mass update and axion forces |
| [field evidence](../stochastic/results/field/README.md) | Archived transfer and refinement checks; one catalogue, no final abundance claim |
| [Fig. 2 evidence](../stochastic/results/paper_figure2/README.md) | Sharp-wall paper-parameter run; archived N64-to-N128 endpoint energy difference 5.83%, not a continuum result |
| [pockets formation kernel](../pockets/formation_potential.h) and [planar methods](../pockets/PLANAR_FIELDS.md) | Possible reuse of coupled forces, exchange diagnostics and small-dimensional checks; different model |

These findings come from source inspection and recorded reports; old checks
were not rerun during this research setup. Existing uncommitted work elsewhere
in the repository was left untouched.

## What the old model contains and omits

`stochastic` samples bubble births, rejects converted birth locations, and
constructs periodic first arrivals. Some experiments use the Poisson sampler;
others use the pinned paper tracer algorithm. They are distinct methods.
For a common prescribed speed, a point converts when a bubble first reaches it.
The axion then sees either a sharp mass switch or a smoothed arrival front.
Smoothing that front does not introduce a dynamical S field.

This is a controlled phenomenological description. It does not evolve wall
stress, acceleration, scalar radiation, collisions, recoil or a plasma. Overlap
is a union of spheres. Conversion is irreversible by construction, so one
arrival time per point suffices. A dynamical S may oscillate or cross back
through a phase boundary; a single arrival map cannot encode its later state.

## A possible mathematical realization, not a chosen model

Use a canonically normalized real S and a fixed axion decay constant f:

    L = (1/2)(partial S)^2 + (f^2/2)(partial theta)^2 - V(S,theta;T)
    V = U(S,T) + f^2 m_a^2(S,T) [1 - cos(theta)].

Here U has a metastable false minimum and a lower true minimum separated by a
barrier. At a chosen reference temperature, call their locations 0 and v.
Require m_a^2(0)=0 and m_a^2(v)=M_a^2 if matching the paper's endpoint masses.
The values between and beyond these endpoints need a physical specification.

One illustrative effective choice is

    U(S) = mu^2 S^2/2 - A S^3/3 + lambda S^4/4 + C
    m_a^2(S) = M_a^2 (S/v)^2.

For this U, mu^2>0, lambda>0, A>0 and A^2>4 lambda mu^2 give additional
stationary points; one must also require the nonzero minimum to be lower
than U(0). C fixes an energy zero. No parameter values have been selected.
The squared coupling is smooth and nonnegative for real S, but is only a
candidate effective model. It is not uniquely dictated by the paper, and its
large-field stability and UV matching must be checked with the full U.
Other S dependences can produce different dynamics despite identical endpoint
axion masses. A real S is consistent with the explicit field content stated
in Appendix A; we need not evolve the heavy radial part of the PQ field unless
the assumed scale hierarchy fails.

With cosmic time t, comoving position x, H = a_dot/a, and optional cosmic-time
friction coefficient eta, the corresponding candidate equations are

    S_ddot + (3H + eta) S_dot - laplacian(S)/a^2
        + U_S + f^2 (m_a^2)_S [1-cos(theta)] = 0
    theta_ddot + 3H theta_dot - laplacian(theta)/a^2
        + m_a^2(S,T) sin(theta) = 0.

With conformal time d tau = dt/a and Hconf = a'/a:

    S'' + (2 Hconf + a eta) S' - laplacian(S) + a^2 V_S = 0
    theta'' + 2 Hconf theta' - laplacian(theta)
        + a^2 m_a^2(S,T) sin(theta) = 0.

The f^2 factor in the S force matters. Rescaling S to S/v also changes its
kinetic normalization; do not drop v^2. Here f is constant: integrating out a
radial field with appreciable S dependence could require a different kinetic
theory. These were the research-stage equations. Their eta=0,
temperature-independent version is now implemented for the selected polynomial
potential and prescribed radiation background; see COSMOLOGY_INPUTS.md.

## Choices that cannot be hidden in the numerics

1. **Nucleation.** A perfectly homogeneous metastable classical field stays
   there without a seed or fluctuations. For a first controlled comparison,
   prescribe birth events and insert finite-radius bubble field profiles with
   declared S and S_dot. This evolves post-nucleation dynamics, not tunnelling.
   A thermal treatment could later compute bounce profiles and a rate such as
   Gamma(T) ~ prefactor exp[-S3(T)/T]. Quantum tunnelling is a different regime.
   Thermal noise would need a specified bath, fluctuation-dissipation relation,
   regulator and convergence treatment; arbitrary grid noise is not a substitute.
2. **Wall motion.** A scalar-only biased potential generally accelerates walls.
   It does not automatically reproduce a constant subluminal speed. Matching
   v_w=0.15 requires a justified friction/transport prescription or an explicitly
   approximate calibrated damping model. Scalar damping alone does not simulate
   plasma hydrodynamics. Collisions and axion pressure can change the speed.
3. **Feedback and nucleation environment.** For nonzero theta, the added axion
   energy changes the full potential, barrier and pressure. The vacua of U alone
   need not be stationary points of V. We must check false-vacuum metastability,
   true-vacuum shifts and whether a theta-independent birth rate is justified.
   Small coupling is not enough without comparing forces and energy scales.
4. **Potential matching.** An S-dependent term independent of theta exerts a
   force on S. Rewriting -chi(S) cos(theta) as chi(S)[1-cos(theta)] requires
   moving -chi(S) into U. It cannot be discarded as a constant. Appendix A
   must be matched at the level of the full potential before assigning forces.
5. **Expansion and thermal bath.** A prescribed radiation background is a useful
   first approximation only if the transition and fields do not dominate the
   energy density. If U depends on prescribed T(t), it exchanges energy with an
   external bath. Thermal free energy is not automatically the bath's complete
   internal energy; a later reheating calculation needs a thermodynamic model.
6. **Resolution.** We must resolve wall thickness, Lorentz contraction, scalar
   oscillations, axion wavelengths and bubble separation in the same box.
   The scalar can be much heavier than the axion, making this harder than the
   current single-field problem. A real S has finite critical bubble radius;
   zero-radius births are an approximation requiring a scale hierarchy.

For periodic boundaries, time-independent V and eta=0, the closed two-field
energy obeys rho_dot + 3H(rho+p)=0 after spatial averaging. In Minkowski space
this becomes total energy conservation. With this potential split, the local
S-induced axion energy source is

    Q_a = f^2 (m_a^2)_S S_dot [1-cos(theta)].

The S sector receives the opposite exchange. Explicit temperature dependence,
friction and inserted bubbles add further sources/sinks; insertion is an
externally prescribed tunnelling event and must be logged. A damped model
requires bath accounting rather than a false closed-system conservation test.

## What this could give us

The first scientific payoff is a validity test of the thin-wall approximation:
does a resolved, evolving mass source reproduce the old axion result when the
wall is thin and feedback is weak? Agreement would strengthen the old model.
A difference must survive numerical controls before being interpreted as physics.

Once that limit is understood, we can ask whether finite wall crossing time,
S oscillations after collisions, scalar radiation or axion recoil change the
axion energy, spectrum and eventual abundance. We can track where axion energy
comes from instead of treating the mass switch as an external source. Useful
ratios include M_a * wall thickness, scalar curvature mass / M_a, H/M_a,
beta/H, critical radius / bubble separation, and axion energy / driving vacuum
energy. For a moving wall, the crossing-time ratio also depends on wall speed
and Lorentz contraction.

The full axion potential energy scale is of order f^2 M_a^2; compare it with
the actual transition pressure and barrier for the selected theta. Feedback
may reduce the driving pressure substantially. This is a reason to audit the
hierarchy before a 3D run, not a prediction that it will occur.

Gravitational-wave production becomes a possible later application through
anisotropic stress, but requires its own tensor evolution and adequate source
physics. S evolution alone does not establish a gravitational-wave prediction,
plasma sound-wave signal or improved axion abundance.

## Proposed order after the research discussion

This was the prospective research sequence. The user has since authorized
the first experiment and its prerequisite checks; later stages remain proposals.

1. Choose the question: establish the old model's domain of validity, or study
   a specified microscopic model. My recommendation is the former first.
2. Write a model contract: U, m_a^2, kinetic normalization, scales, bath,
   expansion, feedback treatment and nucleation/insertion prescription. Check
   full-potential vacua and accessible hierarchies analytically before coding.
3. Later validate homogeneous coupled evolution and energy exchange against an
   independent ODE. Then test a resolved planar wall and a single expanding
   bubble, separating time, space and boundary refinements.
4. Later compare a matched prescribed-front axion calculation with dynamical S
   for one wall/bubble, first in a spectator limit, then with feedback.
5. Only then attempt collisions and small stochastic ensembles. Reuse birth
   events for paired comparisons where physically justified, but do not claim
   identical arrival maps or false-vacuum exclusion once growth differs.
6. Relic-density or spectrum statements need late-time diagnostics, ensemble
   errors, finite-volume and continuum controls. They are not pilot outputs.

For comparison to `stochastic`, preserve c=1, comoving lengths and explicit
cosmic/conformal labels. Its radiation normalization is a=H_* tau and
t=H_* tau^2/2. Its physical rate maps to q=Gamma*a^4; explicit false-vacuum
rejection must not be combined with a second false-fraction factor. Its
constant-speed phantom-equivalence argument does not automatically extend to
finite-radius, accelerating, colliding or axion-dependent bubbles.

## Questions for the next research conversation

- Does the supervisor want an effective real order parameter or a specific
  finite-temperature model based on Appendix A?
- Is the immediate aim to test thin-wall accuracy, or to measure feedback and
  collision effects? Those require different parameter hierarchies.
- Can the first dynamical-wall calculation use a scalar-only accelerating wall,
  or must it reproduce the paper's slow terminal speed?
- Should nucleation remain phenomenological initially, with bounce-derived
  nucleation treated as a later project?

Suggested first target, pending discussion: one seeded bubble in a deliberately
accessible effective model, with a spectator axion control. This isolates what
evolving S adds before combining nucleation, plasma physics and 3D statistics.
