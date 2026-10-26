# Field-only balanced-energy pilot, 2026-10-17

Contract: [FIELD_GOAL.md](../FIELD_GOAL.md), recorded before launch.
Driver: [run_balanced_pilot.py](../run_balanced_pilot.py).
Raw provenance, initial/final states, slices and histories are in the fresh
ignored directory `outputs/pockets/balanced-pilot-01/`.

This is the unchanged compiled phi^2 coupling model, with the gas energy
raised from 1554.9 to 2680.8 using the initial thin-wall pressure estimate.
No particles, thermalisation, expansion or new mass function were introduced.
Time is static cosmic/conformal time (a=1); all quantities are dimensionless.
N=128, L=64, dx=0.5, dt=0.0125, R0=20, w=2, lambda=2,
DeltaV=0.01, m_TV=2; bound modes 0.15<=omega<=0.6, seed=0.

## Outcome

The solver completes t=80. R decreases to 16.687 (also its minimum),
compared with the original low-pressure trajectory's R near 11 by t=100.
These are different initial energies; this is not a controlled lifetime
comparison and does not establish stability.

- Maximum relative energy drift: 1.011e-5, passes the 5e-4 pilot gate.
- Initial-state transfer: exact.
- Stored far axion energy at t=80: 3.176% of initial axion energy.
- Minimum measurable central-row wall width: 4.208 cells.
- Wall width unmeasurable by the two-crossing diagnostic at t=24,60,68,72,76.

**Numerical pilot gate fails.** The unmeasurable slices were counted as
failures, not omitted to obtain a passing width criterion. Multiple crossings
can result from wall deformations or dilaton radiation; their cause has not
been separated. The box also permits periodic radiation to return within
this duration. The far energy is storage, not an irreversible escape fraction.

No relic pocket, isolated lifetime or mass-radius scaling has been established.

## Next work

Inspect the wall crossings and replace the single central-row diagnostic
with a geometrical wall-resolution measure that records extra crossings and
cannot silently discard failures. Refine time and space independently for
this fixed initial continuum state. Check box dependence at fixed dx and
physical parameters before interpreting long-time loss. Consider a different
mass function only as a separately documented field model comparison.

The unchanged compiled backend passes all 25 checks in
`outputs/pockets/validation-goal-01/`. Eight pocket reference unit tests pass;
the optional backend unittest is skipped because its environment variable is
unset, but the backend suite was run directly.
