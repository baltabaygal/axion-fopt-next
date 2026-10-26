# Controlled finite-time planar wall reversal

2026-10-19. The coherent axion field reverses the initially contracting
planar wall, and the response persists under separate temporal/spatial
refinement and displacement of the second wall. This completes the listed
controls for the finite-time pilot in PLANAR_PRESSURE_PILOT.md. It does
not demonstrate a long-lived pocket, equilibrium, physical mass/radius
scales, or cosmological evolution.

The pilot's canonical fields, spectrum, dimensionless scales, Minkowski
time, energy normalization and reversal gates are unchanged. Each control
uses a matched empty-wall evolution. Refinements hold physical initial data
fixed, with no spectrum or velocity fitting. All native runs use one thread.

| Control | Result | Declared requirement |
|---|---:|---:|
| Time field-difference ratio, dt=.025/.0125/.00625 | 4.00490 | >=3 |
| Space field-difference ratio, N=2048/4096/8192 | 15.32565 | >=12 |
| Finest spatial trajectory difference from baseline | .000531672 | <.02 |
| Sampled turning-time change | 0 | <.5 |
| Displaced-wall local scaled field difference | .00000714694 | <.003 |
| Displaced-wall trajectory difference | 2.84e-13 | <.002 |
| Displaced-wall axion-energy-ratio difference | 2.30e-11 | <.0002 |

Every refinement retains the sampled turn at t=102.75 and positive contour
velocity through the final 76.5 units. The second-wall test extends the
physical domain from [0,512] to [-256,512], moving that wall from 0 to -256
while preserving the initial target/packet core exactly. Its prospectively
chosen final comparison core [64,288] excludes the intentionally changed
second-wall region. This is one finite-time boundary sensitivity test.

Matched energy closure and exchange accounting pass their unchanged .002
gates. Refinement runs cost 14.14 solver seconds, and the displaced-wall
pair costs 4.24 seconds. No new 3D evolution was launched.

Evidence lives in fresh directories under outputs/pockets/:
planar-pressure-refinement-01 and planar-pressure-second-wall-01, each with
contracts, histories, final fields and source/build provenance. The compact
index planar-pressure-controlled-01 hashes those reports and the original
pilot report. Runners are refine_planar_pressure.py and
run_pressure_boundary_control.py. The original pilot's pending-validation
status remains historical; its raw report is not rewritten.

Next: prepare a positive-vacuum-bias contraction experiment with distributed
trapped field and a matched empty control. The present degenerate 1D walls
have neither a volume bias nor spherical curvature, so their rebound alone
cannot establish pressure-supported pocket equilibrium. Survival, scaling
and FRW expansion remain open parts of FIELD_GOAL_REVISED.md.
