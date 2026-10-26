# L20 transient-conversion refinement

Declared after the L20 exploratory run and before these refinements. Its
all-cell S>=.5 interval was [28.16,29.24], lasting 1.08 time units, with
recrossing afterward. It did not meet the original sustained >=2 criterion.
The existing report remains goal_reached=false. Energy/statistical gates
passed, but its pre-collision spherical-width gate could not apply: the
initial .1 contour already lies outside L/2, so no isolated spherical wall
window exists at L20. This is an initial periodic-interaction configuration.

Do not erase that failed/inapplicable gate or call the original run a complete
validated success. Test the different, narrower claim of a numerically stable
transient all-cell interval directly using two controls at L20, unchanged
physics and continuum seed, global history every .04, slices .2:

- Spatial: N96, dx=20/96, dt=.01, t=0..30.
- Temporal: N80, dx=.25, dt=.005, t=0..30.

The original baseline N80,dt=.01 is compared only through t30. One CPU
thread per subprocess, no simultaneous field runs, timeout 180 seconds each.
Keep full initial/final states and phase/energy diagnostics. Same reference
validation and passive global-phase checks; source/build hashes preserved.

Gates: finite fields, exact initial transfer, independent phase statistics
(fraction exactly, moments 1e-12), independent endpoint energies 1e-11;
energy drift / initial FV excess <.005 and exchange residual <.01. No
pre-collision radius or thickness is asserted in this deliberately interacting
periodic state. Completion means precisely the same global S>=.5 threshold.
Both refined runs must show an all-cell interval; start/end shifts <=.4 for
space and <=.08 for time versus the baseline. Axion energy at t30 changes
by <5% spatially and <1% temporally. Retain the outcome if any gate fails;
two grids and two timesteps are not a continuum proof.

These controls assess transient conversion only. They do not make the >=2
sustained goal pass, imply relaxation to S=1, or authorize damping or a scan.
