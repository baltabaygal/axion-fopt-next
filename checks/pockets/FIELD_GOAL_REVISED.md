# Accepted field-only formation programme

Accepted by the user on 2026-10-18. This supersedes the static-pocket-first
scientific sequence in FIELD_GOAL.md. Previous results and failed gates
remain archived evidence, not a reproduction of the formation paper.

Using coupled axion–dilaton fields only, test whether the formation and
survival mechanism of arXiv:2407.01676 survives field dynamics. First
establish an accessible, controlled parameter hierarchy; then study
contracting-pocket compression, rebound, leakage and damping. Test
equilibrium pressure and mass–radius relations where applicable. Finally
include radiation-dominated expansion and bubble-collision formation,
distinguishing measured results from assumed cosmological inputs.

## Ordered work and evidence

1. Model/hierarchy audit. Specify canonical kinetic terms, vacuum potential,
   exponential axion mass and periodic axion potential, using explicit
   controlled approximations. Validate derivatives and independent evolution.
   Report wall-frame reflection, post-reflection frequencies, wall/penetration
   resolution, initial stress anisotropy, occupations and computational cost.
   Do not claim an inaccessible Fig. 4 hierarchy from a different regime.
2. Contracting-pocket experiment. Use dynamical walls and causal field
   propagation. Measure R(t), energy/stress profiles, compression, wall
   reversal, escape and exchange against matched controls. Refine time,
   space and boundaries separately. A reduced-dimensional field reference
   must be labelled and cannot establish 3D stability or formation efficiency.
3. Survival and equilibration. Measure damping and net energy transport;
   do not require damping to reproduce a paper result whose damping time
   was not established. Retain t>=20 R0, perturbation and numerical-control
   requirements before using the label long-lived. Oscillating survival,
   equilibrium, thermalisation and numerical transients are distinct outcomes.
4. Pressure and mass–radius tests. Compare stresses and full energy budget
   across several initial radii/energies. Test isotropy and relativistic
   equation of state instead of imposing them during compression. Audit
   finite-tension algebra before using it as a quantitative gate.
5. Expansion and formation. Add prescribed radiation domination with a
   continuity check, then bubble collisions, initial velocities, formation
   efficiency and population evolution. A stationary single pocket cannot
   determine the cosmological alpha, epsilon or abundance.

No particles or N-body implementation. Use the actual compiled field
backend for lattice evidence; retain independent references and source/build
provenance. Declare gates before physical comparisons. Negative or
inconclusive outcomes remain valid findings, not permission to relax gates.
Physical masses/radii require a specified scale map and formation inputs.

Immediate work: [paper/model audit](PAPER_FORMATION_AUDIT.md). Existing
memory improvements may support later runs but do not dictate parameter
choices. The completed toy pair and all its failures remain unchanged.


## Current compute constraint and checkpoint

The user requested light computation on 2026-10-18 because the laptop is
running other work. Continue source work, analytic references and tiny tests;
defer production simulations and substantial builds until resources are
available. The standalone C++ exponential/cosine kernel now agrees with the
independent reference. See [integration checkpoint](FORMATION_INTEGRATION.md)
for completed checks and required backend/diagnostic integration. This does
not yet constitute a compiled field evolution of the new model.

A standalone homogeneous C++ evolution now passes independent DOP853,
timestep-refinement, total-energy and internal-exchange checks over t=0..4
(`outputs/pockets/formation-homogeneous-01`). This prepares the reference for
the future framework comparison; it is not a spatial field-backend result.

Optional formation selection is now integrated in CPU backend source, full
mixed curvature diagnostics, and both energy monitors. A framework syntax
check and independent fixed/moving continuity tests pass. Rebuilt-framework
runtime validation and old-model compiled regression checks remain pending;
no production runs were launched under the current compute constraint.

Compiled fixed/moving monitor agreement and instantaneous continuity now pass
for formation and two polynomial variants on a 12^3 state, using the actual
C++ monitor headers (`outputs/pockets/formation-monitor-agreement-02`). The
shared scalar dispatch also passes independent-reference tests. Full framework
evolution, integrated accounting and archived-run regression remain pending.

Actual framework evolution now passes a homogeneous 4^3 comparison with
independent DOP853, PV2 timestep refinement, energy, mixed-Hessian frequency,
uniformity and parameter-rejection checks. A serial O0 driver build reused
existing libraries; archived polynomial runs match exactly at N=8 for p=2/4,
stiffening 0/12. Evidence is in formation-backend-homogeneous-01 and
formation-polynomial-regression-01 under outputs/pockets. These tiny checks
respect the laptop constraint; physical pocket experiments remain deferred.
Next: independent inhomogeneous dynamics and integrated monitor continuity.

Small inhomogeneous actual-framework evolution now agrees with independent
fixed-stencil DOP853 and shows second-order temporal refinement. Integrated
fixed/moving sector continuity, independent output-cadence auditing and
periodic translation controls pass (`formation-backend-inhomogeneous-01`).
N=8 is a diagnostic state only, not resolved pocket physics. Next choose
accessible wall-frame reflection and coupled-field initial-data hierarchies;
production simulations remain deferred under the laptop constraint.

Analytic candidate audit now records initial moving-wall resolution,
wall-frame reflection and repeated-encounter failure, finite-mode counts,
finite-amplitude dilaton backreaction and isolation limits. See
[prospective hierarchy](FORMATION_CANDIDATES.md) and
`outputs/pockets/formation-candidate-audit-04`. No candidate grid was run.
Local scattering and initial-profile preparation precede production runs.

An independent linear planar scattering ODE now passes analytic constant/step
controls, endpoint variation and tolerance refinement. It quantifies strong
transmission at the candidate upper-edge second encounter for an unchanging
wall speed; see FORMATION_SCATTERING.md and formation-scattering-reference-01.
This zero-bias weak-field reference is not a coupled pocket survival result.

A native NDim=1 CosmoLattice coupled-field driver now passes independent
inhomogeneous ODE, PV2 refinement, energy/exchange and periodic translation
checks (formation-line-validation-01). It shares the production scalar
potential, preserves the 3D driver, and avoids cube allocation for local
scattering. See PLANAR_FIELDS.md. This is a planar reference capability and
cannot replace 3D survival/scaling tests. No production run was launched.

Deterministic positive-k weak-packet and periodic kink-pair data now pass
initial directionality/localization/boost audits, with fixed continuum
coefficients on two grids and matched empty-wall controls. Eight input
cases are prepared but not run (planar-scattering-preparation-02). The
continuum kink's lattice residual is measured, not assumed absent. Next:
empty-wall/runtime control and validated local scattering diagnostics.

Short t<=4 native planar empty-wall controls pass on two grids, with improved
continuum agreement and exact zero axion leakage (planar-empty-control-01).
Independent wall geometry and moving-window sector-current diagnostics pass
small tests. These controls cost 2.30 seconds and do not certify the full
packet window. Next expose local diagnostics in the compiled line output.

Compiled planar local window energies, currents, exchange and advection now
pass independent endpoint agreement, integrated continuity/time refinement,
passivity and periodic translation checks (planar-window-validation-01).
The new driver remains one-thread and a separate binary; no full scattering
run was launched. Next: a bounded weak-packet case with matched controls.

A first single-grid stationary planar packet/control pilot ran in 29.20
one-thread seconds under a 45-second budget. Axion energy largely returns
to the FV side, but the detector endpoint-coverage gate fails (only 93.3%
enclosed). The finite-band continuum reference similarly encloses 93.7%,
so this is a detector/arrival-window issue, not measured 6.7% loss. The
original failed gate is retained; no scattering coefficient is validated.
See PLANAR_STATIONARY_PILOT.md. Next redesign detector coverage using the
reference before independently refining time/space/amplitude/boundaries.

The detector geometry has been redesigned prospectively using the finite-band
reference, preserving the 98% coverage requirement. Reference endpoint
coverage is 99.745%; offline old-endpoint coverage is 99.717%, without
altering the original failed report or claiming a new flux measurement.
The bounded runner accepts the hashed design for a future run; no solver
ran during this step. See PLANAR_DETECTOR_DESIGN.md.

The revised-detector actual stationary pilot now passes all predeclared
accounting/coverage gates (planar-stationary-pilot-02); physical scattering
validation remains false pending independent refinements. An O2 build with
no fast-math repeats the pair in .921 single-thread seconds, with exact
endpoint fields and axion histories. Full-history bitwise equivalence fails
on a single 1.11e-16 dilaton/total diagnostic entry and is retained as failed,
with a separate precision audit. See PLANAR_PILOT_ACCEPTED.md. The lower
cost enables time/space/amplitude/boundary controls while respecting the
laptop constraint; no 3D production run has started.

Time/space/amplitude stationary planar controls now pass within 24.925
single-thread solver seconds: successive field-difference ratios are 3.999
in time and 3.886 in space; half-amplitude scaled axion waveforms differ by
1.60e-4 and matched recoil scales by 4.000. Apparent tiny TV energy falls
~16x per grid halving, so no physical ppm leakage rate is claimed. Boundary
control and physical-scattering validation remain explicitly false; see
PLANAR_REFINEMENT.md. Next displace the second wall while preserving the
local target-wall and packet data, then test moving walls.

The finite-time stationary weak-packet planar benchmark now passes its listed
controls, including a displaced-second-wall pair (7.34 one-thread seconds)
and a direct continuum waveform test improving ~4x per grid halving to
0.429% of fixed incident scales. Machine index: planar-stationary-benchmark-01.
This is a bounded local reference result, not a precise asymptotic loss rate
or 3D survival. See PLANAR_STATIONARY_BENCHMARK.md. Next: moving-wall Doppler
and energy-exchange benchmarks before coupled compression experiments.

A first native moving-wall pilot passes declared energy/Doppler/accounting
checks in 4.35 one-thread seconds: measured lab energy gain 1.84225 vs
1.84392 reference, energy-weighted frequency .93559 vs .93603, with matched
dilaton energy supplying the gain. The corrected moving reference passes
initial-packet purity after spectral-alias detection; prior failures remain
recorded. Waveform phase error is still .225 in initial scales, so moving
scattering is explicitly unvalidated pending independent refinements. See
PLANAR_MOVING_PILOT.md. Next: moving time/space/amplitude/boundary controls.

Moving time/space differences now converge by factors 4.000/3.994 and the
half-amplitude control passes, but the direct 1% waveform gate still fails
at 5.568% in fixed incident scales (planar-moving-refinement-01). Empty-wall
position drift scales quadratically and explains much of the phase error;
its diagnostic effective-speed reference does not replace the failed gate.
A 6.52-second moving second-wall pair has negligible energy/recoil impact,
but its broad final region incorrectly includes the deliberately swept
second wall and fails. A separate local-interior audit is small; future
runner contracts now distinguish initial overlap from final interior.
Both original failures are retained. See PLANAR_MOVING_REFINEMENT.md.
No large 3D run was launched; next improve moving-wall spatial accuracy
within the compute constraint, without fitting away the fixed-speed error.

A fourth-order planar spatial reference now passes force/Hamiltonian,
long-bond moving-window continuity and Fourier-order checks. The sampled
moving-kink residual is ~35x smaller at fixed dx=.125; no compiled evolution
or simulation changed. See PLANAR_SPATIAL_ACCURACY.md and
planar-fourth-order-reference-01. Next integrate an optional owned adapter
hook with matching energy/current/frequency bound and two ghost cells,
then retest the same physical moving benchmark and unchanged 1% gate.

The optional order-4 planar operator is now compiled with two ghost cells,
matching signed-bond energy/current and frequency guard, via an owned overlay
that leaves vendor files untouched. Independent small backend/window checks,
actual fourth-order spatial convergence and default order-2 exact regression
pass. The same N=4096,dt=.0125 moving packet now passes the unchanged fixed-u
1% waveform gate at .108%, versus 22.5% for order 2, in 5.03 one-thread
solver seconds. See PLANAR_ORDER4_INTEGRATION.md. This first new-scheme pilot
is not full moving validation: next repeat physical time/space/amplitude and
corrected-region boundary controls before contracting-pocket studies.

The fourth-order finite-time moving-wall benchmark now passes the complete
listed local controls: temporal ratio 4.000, spatial ratio 15.538, finest
fixed-reference waveform error .0862%, weak-amplitude accounting, and a
prospectively corrected-region second-wall comparison. New controls cost
27.07+6.72 single-thread solver seconds under separate short caps; no 3D
run. Evidence index: planar-moving-benchmark-01; summary:
PLANAR_MOVING_BENCHMARK.md. Old second-order and bad-region failures remain.
Next test finite-amplitude coupled-field pressure slowing/reversal, without
assuming a stable pocket or replacing the still-unmet 3D/scaling/FRW goals.

A finite-amplitude native planar pressure pilot observes sustained wall
reversal: A=.16 changes final contour velocity from ~-.299 to +.235 and
advances 4.71 units in the final 20; A=.08 slows but does not reverse.
Three one-thread runs cost 6.55 seconds and pass declared energy, exchange,
geometry and resolution checks. The pilot has degenerate vacua, no stable
pocket claim, and turnaround validation remains false pending refinements
and displaced-wall controls. See PLANAR_PRESSURE_PILOT.md and the inspected
trajectory plot in planar-pressure-plots-02. Next validate this nonlinear
response, then bias-driven/distributed-field contraction before 3D survival.

The finite-amplitude planar reversal now passes separate time and space
refinements (ratios 4.005 and 15.326) and a displaced-second-wall check.
All four refinement pairs plus the boundary pair cost 18.38 one-thread
solver seconds. Turning time remains t=102.75; the finest spatial trajectory
differs by at most .000532. See PLANAR_PRESSURE_CONTROLLED.md and evidence
index planar-pressure-controlled-01. This establishes a controlled finite-time
degenerate-planar response, not survival or equilibrium. Next prepare
positive-bias/distributed-field contraction with matched empty controls and
prospective energy/turnaround/leakage gates; retain the laptop compute limit.

Positive-bias initial data are now explicitly separated from exact degenerate
kink data. Independent force projections at three resolutions recover
F=-DeltaV to <1e-12 relative for DeltaV=.0002. Full packet energy .09146
is comparable to initial false-volume bias energy .0512 per area. Three
lightweight tests pass; no biased solver run was launched. See
PLANAR_BIASED_PREPARATION.md and planar-biased-seed-audit-02 (audit-01
omitted axion gradient energy and is superseded). The next bounded pilot
uses matched empty/packet data and unchanged numerical/reversal gates;
distributed trapped-field and full survival/scaling/FRW goals remain open.

The positive-bias pair now runs in 3.99 one-thread solver seconds. At
DeltaV=.0002 the empty wall accelerates inward to v=-.38887; A=.16
reverses at sampled t=104.5 and finishes at v=+.16022. Energy/exchange,
geometry, width and exact-zero gates pass. This is a single-grid biased
planar rebound, not long-lived survival. See PLANAR_BIASED_PRESSURE_PILOT.md
and planar-biased-pressure-pilot-01. Next refine this response and check
second-wall sensitivity before distributed trapped-field contraction.

The biased rebound passes separate time/space refinement and displacement
of the second wall. Difference ratios are 4.00375 and 15.07441; sampled
turn time remains 104.5 and finest trajectory difference is .000602.
The displaced-wall core [80,288] was declared before evolution, excludes
the biased partner endpoint, and passes with local scaled difference
7.91e-6. Controls cost 20.59+6.39 one-thread solver seconds. Evidence:
planar-biased-pressure-controlled-01; PLANAR_BIASED_PRESSURE_CONTROLLED.md.
Next move from a one-sided impulse to distributed trapped fields and
repeated compression with moving-window energy/flux accounting. Do not
interpret another rebound as equilibrium; full survival/scaling/FRW and
formation goals remain outstanding. No further simulations in this batch.

Distributed standing-field initial data now fill a symmetric biased FV slab,
with explicit continuum modes/phases, taper and canonical momenta. The
linear phase-average pressure scale equals DeltaV; this is not exact wall
balance. Three initial-data tests pass with fourth-order energy convergence.
The tapered initial spectrum has ~9.77e-5 linear energy above the boosted
threshold; no slab evolution has been run. See
PLANAR_DISTRIBUTED_PREPARATION.md and planar-distributed-seed-02. Next
measure both-wall compression and window flux using native fields with
a matched empty control; repeated interactions require crossing times ~260.

The first native distributed-field slab pair now reaches T=320 in 8.37
one-thread solver seconds. Filled length decreases 256->207.95 versus
163.92 empty; final right-wall velocities are -.1067 and -.2353. Both
continue contracting, so no equilibrium is claimed. Declared energy,
exchange, stationary-window, endpoint geometry/symmetry and resolution
gates pass. See PLANAR_DISTRIBUTED_PILOT.md and planar-distributed-pilot-01.
Report serialization was repaired and analysis recovered from saved fields
without rerunning either solver. Next refine and check exterior-buffer
sensitivity, then moving-boundary retention and longer interactions.

The distributed slab's slowing persists at coarse time and coarse space
levels with unchanged numerical gates. Filled endpoint length differences
are .000258 (dt doubled) and .009717 (dx doubled). New pairs cost 7.84
one-thread solver seconds; two-level comparisons do not prove order.
See PLANAR_DISTRIBUTED_PILOT.md and planar-distributed-coarse-comparison-01.
The runner now accepts n/dt with fixed cadence and physical seed. Fine
levels and exterior-buffer/moving-retention controls remain open.

A saved-endpoint wall-following audit exposes ~3.12% of current axion energy
outside a window around the actual slab interfaces (padding=4), despite
tiny flux through the larger fixed window. Coarse-space/time values are
3.112%/3.120%; padding 2..8 gives 3.06..3.54%. This includes penetration
and transition-region energy, not proven irreversible escape. No solver
run was used; two independent window-derivative tests pass. See
PLANAR_SLAB_RETENTION.md and planar-distributed-endpoint-retention-01.
Next obtain time-resolved moving-window exchange/advection and buffer
sensitivity before long-time survival; fixed-window retention is inadequate.
