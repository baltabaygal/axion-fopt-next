# Moving-wall refinement: convergent, with a retained waveform failure

2026-10-19. The actual NDim=1 moving-wall time/space/amplitude suite completed
in 25.10 one-thread solver seconds, within its 35-second wall-time budget.
Evidence: `outputs/pockets/planar-moving-refinement-01/report.json`, declared
contract, per-case raw fields/histories and executable/source provenance.
No 3D simulation ran. This report records incomplete validation, not a
promotion of the pilot to a fully controlled moving benchmark.

Fixed physical parameters: L=512, t=0..180, initial approach speed .3,
degenerate lambda=8, mFV=.1, mTV=1, fa=10, same continuum packet coefficients
and detector as PLANAR_MOVING_PILOT.md. No per-grid speed or phase was fitted.

## Passing controls and unresolved gate

At N=4096, dt=.025/.0125/.00625, successive full-field differences are
.00332193/.00083042 in fixed component scales (1,.001,1,.0005), giving a
ratio of 4.0003. The finest temporal difference passes the .002 gate.

At dt=.0125, N=2048/4096/8192, successive differences are .676543/.169378,
giving a ratio of 3.9943. Fine samples are shifted to the actual coarse
cell centres with the tested periodic Fourier interpolation. Halving the
amplitude at N=4096 changes the rescaled axion waveform by .00040191,
below the .003 gate. All individual energy/frequency/accounting/coverage
pilot checks pass, including exactly zero axion controls.

| N at dt=.0125 | Lab axion energy gain | Energy-weighted frequency | Direct waveform error |
|---:|---:|---:|---:|
| 2048 | 1.837314 | .934290 | .89947 |
| 4096 | 1.842246 | .935587 | .22506 |
| 8192 | 1.843457 | .935894 | .05568 |
| Continuum weak-field reference | 1.843921 | .936026 | — |

Direct errors use fixed incident field/momentum scales, not a fitted packet
phase or the amplified final peak. They decrease approximately fourfold,
but the declared finest direct-error gate is .01. It FAILS at .05568.
The suite exits nonzero, retains `declared_checks_passed: false` and
`moving_scattering_validated: false`. Integrated energy agreement alone
does not override the failed waveform check.

## Empty-wall phase diagnosis, without changing the reference gate

`audit_moving_phase_error.py` reads saved empty-wall endpoints; it launches
no solver. Evidence: `outputs/pockets/planar-moving-phase-audit-01`.
For N=2048/4096/8192, the empty wall's endpoint displacement relative to the
continuum x=256-.3t is .379606/.093199/.023094. Its independently inferred
average approach speed tends to .3 under refinement.

A counterfactual weak-field reference using that average endpoint speed,
with interpolated steady scattering amplitudes, reduces the corresponding
waveform errors to .34717/.08716/.02147: a factor ~2.59 in each case.
This supports wall drift as a major phase-error contribution. It is not
a measured full trajectory, and dispersion/transients remain. Crucially,
it does not replace the original fixed-speed reference or make its .01
gate pass. Further spatial accuracy or refinement remains needed.

## Moving second-wall test and comparison-region issue

A separate packet/empty pair displaced the second wall from x=0 to -256
in physical coordinates, keeping the target x=256, packet, detector,
dx=.125, dt=.0125 and t=180 fixed. It took 6.52 seconds on one thread,
below its 15-second budget. Initial fields and momenta match exactly in
physical [32,288]; far-tail taper changes only 1.10e-17 of incident energy.
Evidence: `outputs/pockets/planar-moving-second-wall-01`.

Energy gain changes by -4.73e-9; matched target recoil changes by only
5.68e-14. Local continuity, detector coverage and the zero control pass.
However, the broad final full-state comparison in [32,288] FAILS, with
scaled difference 1.000004 near x=37.19 in phi. That region is swept by the
original second wall, which ends at x=54, while the displaced wall ends
near -202. The two experiments intentionally have different phases there.
The original boundary report retains its failed check and is not called a
passing full comparison.

A read-only region audit excludes eight nominal widths from both known
nominal wall paths, giving final local interior [62,194]. This region is
derived from initial geometry, speed and end time, not by fitting the
packet output. Its full-state difference is .00013974, below the unchanged
.003 local tolerance. Evidence:
`outputs/pockets/planar-moving-boundary-region-audit-01`. This is a scoped
endpoint audit, not a retrospective replacement for the failed broad gate.
The runner now declares separate initial-overlap and final-interior regions
for future execution; unit tests check the geometry. No duplicate run was
launched merely to change the report's colour.

The moving benchmark therefore remains incomplete. The important next gate
is reducing fixed-speed waveform error without fitting away wall drift,
while retaining the corrected local boundary scope prospectively. No
three-dimensional pocket lifetime, scaling, expansion or formation result
follows from these planar checks.
