# Longer single bubble: declared contract

2026-10-19. The user requested longer single-bubble evolution and all input
parameters. Retain the old t<=6 evidence. Target t=18, with the same action,
potential, units and seed as EXPERIMENT.md. No expansion, friction, bath,
stochastic births or imposed speed. Full reciprocal cosine coupling remains.

Physical inputs: c=1, canonical S, canonical A=f theta, excursion v=1;
lambda=2, DeltaV=.06, M=.5, f=.2. S0=(1-tanh((r-8)/2))/2 at the box centre;
S_dot=0, theta0=1, theta_dot=0. Empty control uses theta0=0. Periodic cubic
box, no random seed. This is an imposed seed, not bounce nucleation. Time
is Minkowski cosmic time and equals conformal time because scale factor=1.

Predeclared cases (each starts from the same continuum initial profile):

| Case | N | L | dx | dt | End |
|---|---:|---:|---:|---:|---:|
| coarse | 96 | 48 | .5 | .04 | 18 |
| main | 128 | 48 | .375 | .04 | 18 |
| time-half | 128 | 48 | .375 | .02 | 18 |
| empty | 128 | 48 | .375 | .04 | 18 |
| larger-box | 160 | 60 | .375 | .04 | 18 |

PV2 integration, nearest-neighbor periodic Laplacian and compatible forward
gradient energy, double precision, one process and one CPU thread. Scalar
history every .08; central 2D slice every .4; full 3D initial/final states.
Fresh run directories, source/build hashes and parameters preserved. A finite
180-second subprocess guard permits longer runs; tool calls yield while they
run so progress can be communicated. Raise only the owned backend grid guard
from N<=128 to N<=192 to permit the boundary control; the field equations and
integrator do not change. Rebuild and repeat the independent small-grid
reference validation before the longer runs.

Prospective gates:

- Reference validation passes on the rebuilt executable, including second-order
  time refinement and snapshot/independent-energy agreement.
- Every initial state transfers exactly; saved states remain finite.
- Energy drift / absolute initial excess above uniform FV energy <.005, and
  independent final snapshot energies / initial total agree within 1e-11.
- Axion integrated-exchange residual / initial axion energy <.01. Empty axion
  remains exactly zero.
- Halving dt changes final radius <.5%, axion energy <1%, and fixed-scale
  [1,f,1,f] full-state RMS <.01. Temporal order was already established for the
  integrator; this longer check measures accumulated error at fixed dx.
- N96-to-N128 endpoint radius changes <5% and axion energy <10%; these are
  two-grid sensitivity checks, not a continuum-order claim.
- Main and empty walls have unique outward .1/.5/.9 central-axis crossings
  and >=5 cells between .9 and .1 at every saved time. Record all crossing
  counts and the first failure time; never silently choose a convenient root.
- L48/N128 versus L60/N160: final local fixed-scale RMS in r<20 <.001 and
  radius difference <.5%. Local comparisons use aligned cell centres. Report
  global axion-energy sensitivity too; a local gate does not certify all cells.

If the wall becomes unresolved, retain the failed run and restrict claims
to its measured resolved window. A predefined refinement option is N160,L48
(dx=.3), dt=.02 plus N192,L57.6 (same dx), dt=.02, with all physical seed
parameters unchanged. Its boundary comparison also uses r<20. Do not reduce
resolution thresholds to obtain a pass. No later automatic extension beyond
t=18 or stochastic/abundance production is part of this request.
