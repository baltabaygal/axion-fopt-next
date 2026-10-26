# Phi^2 / phi^4 field coupling pilot

Contract: [COUPLING_EXPERIMENT.md](../COUPLING_EXPERIMENT.md), declared
before implementation or runs. Raw runs: `outputs/pockets/coupling-pilot-01`.
The prior phi^2 source and executable were preserved under
`outputs/pockets/phi2-source-baseline-01` before changing the implementation.

Both coupled real-field runs complete t=80, N=128, L=64, dx=.5,
dt=.0125; lambda=2, DeltaV=.01, m_TV=2, R0=20, w=2, a=1.
Each has its own .15<=omega<=.6 bound-mode Gaussian high-occupation
initial field, seed=0, initial lattice axion energy 2680.825731.

The new potential is U(phi)+m_TV^2 phi^p a^2/2, p=2 or 4. It is an
explicit toy mass-function family, not an implementation of the paper's
instanton mass. In p=4 the coupling's phi curvature vanishes at phi=0.
The turning surfaces and mode functions also change.

| Measure at t=80 | p=2 | p=4 |
|---|---:|---:|
| Real harmonic mode count | 93 | 98 |
| Gaussian amplitude | 5.61508 | 5.42801 |
| Final radius | 16.6873 | 18.6006 |
| Radius loss | 3.3130 | 1.3997 |
| Axion-sector depletion | 24.355% | 13.139% |
| Stored far axion energy / initial E_a | 3.176% | 1.883% |
| Relative total-energy drift | 1.01e-5 | 6.38e-6 |

The p=2 trajectory reproduces the original balanced pilot within the
predeclared thresholds (including all shared history samples).
Initial-state transfer, completion and total-energy gates pass for both.

**The improvement criterion fails.** Radius loss improves by more than a
factor two, but depletion improves by only 1.85, below the factor-two gate.
Both also fail the all-slice diagnostic: p=4 has 5 ambiguous contour-edge
brackets in total, although its minimum measured in-plane width is 5.06
cells; p=2 has 22 ambiguous brackets and minimum measured width 3.30 cells.
Missing brackets were not discarded to obtain a pass. Full 3D resolution
and a long lifetime are not established by these screening runs.

The reduced depletion motivates further investigation, but does not show
that removing false-vacuum curvature alone produces a relic. Different
spectra and mode counts confound a one-parameter interpretation. In
particular, gas frequencies remain comparable to the dilaton scale: sums
of two initial axion frequencies can exceed the exterior dilaton mass
sqrt(lambda/2+12 DeltaV)=1.0583. Whether that coupling channel dominates
the observed radiation has not yet been isolated.

## Verification

`validation-coupling-02` passes 40 compiled-model checks, including an
independent DOP853 reference for p=4, second-order time refinement (ratio
4.00004), all energy-region allocations and a finite-difference check of
the coupled-field Hessian frequency guard. The full repository unittest
run passes 58 tests with 11 optional backend tests skipped; the pocket
backend was validated separately rather than through its skipped wrapper.

Next: isolate the axion/dilaton timescale hierarchy with a predeclared
lower-frequency field pilot, record the accompanying state-count change,
then refine a candidate before committing to the ten-crossing-time run.
The goal remains static long survival, scaling, expansion survival, and
formation; no stage is marked achieved by this pilot.
