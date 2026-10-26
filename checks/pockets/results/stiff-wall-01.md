# Stiff-wall paired grids and translating-pocket energy

Contract: [STIFF_WALL_EXPERIMENT.md](../STIFF_WALL_EXPERIMENT.md).
Raw fields, histories, initial spectra and original reports are in
`outputs/pockets/stiff-wall-01`. Both compiled Serial runs complete t=80.
The original executable was kept unchanged until both were terminal.

Static a=1, R0=20,w=1,lambda=8,DeltaV=.01,m_TV=2,phi^4 coupling,
dt=.0125,L=64; N=192 and 256. One continuum seed=0 Gaussian state with
99 real harmonic modes (.1<=omega<=.6), radial h=.01,r_max=31,r_cut=30,
amplitude 6.87900247 fixed on the N=192 grid. Initial gas energy on that
reference grid is 4356.341813. No thermalisation is imposed.

## Histories agree, but the original candidate gates fail

| At t=80 | N=192 | N=256 |
|---|---:|---:|
| Final radius | 18.7696 | 18.7825 |
| Minimum radius | 18.7288 | 18.7362 |
| Axion-sector depletion | 6.728% | 6.777% |
| Far stored axion energy / E_a(0) | 1.643% | 1.706% |
| Global energy drift | 7.84e-6 | 7.82e-6 |

The paired history gates pass: maximum radius difference .01399, maximum
far fraction difference .000646, and depletion difference .000491. The
machine-readable comparison is `stiff-wall-comparison-01`. Two grids are
a consistency check, not a continuum limit or long survival result.

The <=5% axion-depletion screening threshold fails. The wall gate also
fails. Final full-3D diagnostics find minimum measurable widths 3.66 and
4.82cells. On the fine grid, two of 114,952 contour-edge points have multiple
half-height crossings along the normal before reaching the .1 threshold;
they were not omitted to obtain a pass. Saved profiles in
`stiff-wall-ambiguities-01` show the folds directly. Both grids also have
ambiguous in-plane points at t=60. No complete wall-resolution pass is claimed.

## Energy flow is measured, but the fixed sphere does not contain the body

The fixed radius 26 monitor's continuity residuals are <7e-5, and its
twice-coarser cadence comparison is <1.4e-4, within declared 1e-3 gates.
Net outward phi/axion energies at t=80 are 181.21/97.46 (coarse) and
184.79/102.28 (fine). These are signed currents, including any incoming
waves, not irreversible escape.

Crucially, the pocket centroid moves 6.07 units by t=80. The false-vacuum
extent about the box centre grows to 27.42, outside the fixed radius 26.
The fixed sphere's 4.50% energy export therefore includes translating
bound material and must not be called total pocket radiation loss.
See `stiff-wall-extent-01`.

A separate periodic-centroid snapshot audit centres the same radius 26 on
the moving pocket. It encloses all false vacuum, with 2.35 clearance at
t=80. Its stored energy fractions are 1.0000, .9970, .9740 at t=0,40,80.
This is instantaneous storage including any waves still inside, not a
moving-volume continuity calculation or a proven rest mass. The audit is
`stiff-wall-moving-budget-01`; its geometry passes periodic translation tests.
The initial global momentum/energy ratio is .0683; the centroid drift is
consistent with a moving initial state rather than an imposed stationary gas.

Consequently, sector depletion, radiation and translation must remain
separate. A late local-radius plateau alone does not prove a small loss rate.
Moving-volume current or rest-frame initial data and box/isolation checks
are needed before the unchanged ten-crossing-time goal can be established.

## Diagnostic precision correction, without changing the old verdict

Both original reports also fail the fixed 1e-12 independent initial-energy
check because of sequential accumulation roundoff on millions of cells.
Exact state transfer passes. The coarse axion discrepancy is 1.53e-12.

Subsequent builds use Neumaier compensated host diagnostic sums only;
field equations and integration are unchanged. A fresh one-step replay of
the exact N=192 input reduces the maximum initial energy error to 1.43e-16
and the axion error to zero, at the same tolerance. The original failed
reports remain unchanged. See `precision-replay-01` and
[DIAGNOSTIC_PRECISION.md](../DIAGNOSTIC_PRECISION.md).

## Supporting numerical work and next steps

The optional threaded build passes all 51 pocket checks, the 8 main adapter
reference checks and invalid-input checks. Matched Serial/OpenMP final
states and histories are bitwise identical in the benchmark; four threads
are about twice as fast as one for that modest case. A Serial-to-OpenMP
checkpoint continuation also agrees bitwise with the uninterrupted field
trajectory. Neither test is a physical lifetime result.

An exactly matched empty-gas control completes under its separate
[contract](../EMPTY_STIFF_CONTROL.md), with an independent radial reference
and an explicitly limited early-resolution window. All early gates pass:
the empty reference crosses R=15 at t=9.75, when the gas case has R=19.469.
See [the control result](empty-stiff-control-01.md). A new vacuum-tail
potential was prepared in [STRUCTURED_WALL_EXPERIMENT.md](../STRUCTURED_WALL_EXPERIMENT.md)
to test a larger exterior mass gap without drastically increasing tension
or narrowing the core. It has since completed compiled reference checks and
paired field runs; see [the negative result](structured-wall-01.md).
No long-lived relic, mass-radius scaling,
expanding survival or collision-formation result is claimed yet.
