# Revised-detector stationary pilot and lightweight optimized backend

2026-10-18. The prospective detector design was used in a new actual
CosmoLattice NDim=1 packet/matched-empty-wall pair. Evidence:
`outputs/pockets/planar-stationary-pilot-02/report.json`. Its contract names
and hashes the independent detector-design report. The original pilot-01
failed coverage gate remains unchanged.

Physical parameters, initial fields, N=2048, L=512, dx=.25, dt=.025 and
Minkowski time 0..220 are unchanged from pilot-01. Only the passive window
geometry changes to centre=104, inner=88, outer=120. Initial and final field
files for both packet and empty control match pilot-01 exactly. The two
new runs used one thread and took 29.11 seconds total, below the unchanged
45-second solver budget. No 3D production experiment ran.

## Diagnostic gates now pass

All five predeclared pilot diagnostic checks pass: exact-zero control,
axion continuity, independent cadence comparison, outgoing crossing >99%,
and final missing packet energy <2%. Their thresholds were not relaxed.

| Measured quantity, relative to initial axion energy | Value |
|---|---:|
| Outward crossing | 0.99989265 |
| Returned crossing | 0.99706416 |
| Final detector energy | 0.99717243 |
| Final total axion energy | 0.99998781 |
| Maximum integrated axion continuity residual | 4.684e-5 |
| Output-cadence difference | 8.179e-6 |

The matched continuum finite-band reference predicts outward/returned
crossings 0.99991270/0.99736022 and endpoint coverage 0.99745006. The
measured endpoint differs by ~0.000278 in incident-energy units. Residual
packet energy outside the window is still explicitly present; the returned
crossing is not silently promoted to a unit-asymptotic reflection coefficient.

This is an accepted diagnostic pilot, not a converged scattering experiment.
The coarse final waveform error against the continuum reference remains
~6.7% of incident amplitude. Time, space, amplitude and boundary controls
are still required. The report deliberately retains
`physical_scattering_validated: false` despite the passing diagnostic gates.

## Optimized build for lower laptop load

`build_planar_backend.py` compiles only the same planar driver against the
existing CPU libraries and records library, source, compiler/link and binary
identity. A separate O2 executable is saved in `build/pockets-line-opt-01`.
No fast-math option is added, no forces or integration algorithm are changed,
and no archived executable is replaced. The build took 3.80 seconds.

Independent small-window validation passes with exactly matching baseline
fields and existing history values (`planar-window-opt-validation-01`).
Four short stationary/moving wall controls pass
(`planar-empty-opt-control-01`). Repeating the full revised-detector pair
with O2 took 0.921 seconds, versus 29.11 seconds for O0, on these modest
one-thread runs. Evidence: `planar-stationary-opt-pilot-01`. This is a
measured cost reduction for these cases, not a general hardware-scaling claim.

The deliberately strict full long-history byte comparison FAILS and is
retained in `planar-optimization-agreement-01`. Both initial/final field
files match exactly, as does the entire empty-wall history. In the packet
history, one entry in each of total and phi_energy differs by
1.110223e-16. All other entries, including every axion energy, current,
exchange, advection and frequency-bound value, match exactly. The diagnostic
difference is 3.11e-11 of initial axion energy. A separate precision audit
records these facts without changing the failed exact-history gate:
`planar-optimization-precision-audit-01`.

Thus full-history bitwise equivalence is not claimed. The optimized build
passes the independent dynamics/accounting checks and reproduces the saved
endpoint fields and axion observables exactly on the long packet case.
This makes the required refinements feasible at low single-core cost.

Next refine time and space separately at fixed physical data, halve the
packet amplitude, and vary the second-wall/boundary geometry with controlled
initial overlap. Then test moving walls and the predicted frequency shift.
Three-dimensional long-lived survival, mass-radius scaling, expansion and
collision formation remain unachieved.
