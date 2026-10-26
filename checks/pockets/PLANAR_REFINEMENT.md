# Stationary planar time, space and amplitude controls

2026-10-18. Native CosmoLattice NDim=1 coupled fields, one thread. Evidence:
`outputs/pockets/planar-stationary-refinement-01/report.json`, predeclared
contract, raw initial/final fields, histories, logs and source/build identity.
Nine new solver runs cost 24.925 seconds total, below the declared 30-second
budget. No 3D grid was allocated. The earlier optimized baseline is reused
by hash rather than rerun.

Physical data remain L=512, t=0..220, degenerate lambda=8, mFV=.1, mTV=1,
fa=10, initial positive-k packet with central omega=.5, envelope width 16
and amplitude .001. The continuum Fourier coefficients and canonical
momenta are identical across grids; they are not retuned to lattice
dispersion. The detector remains centre=104, inner=88, outer=120, with
output interval .5. All controls have matched exact-zero axion wall runs;
the half-amplitude case reuses the identical N/dt empty control.

## Separate refinements

Time is refined at fixed N=2048: dt=.025,.0125,.00625. Successive final
full-field differences are .00061942 and .00015489 in the fixed component
scales (1,.001,1,.0005), a ratio of 3.9991. This passes the declared
factor-3 order gate and .002 finest-difference gate.

Space is refined at fixed dt=.0125: N=2048,4096,8192 (dx=.25,.125,.0625).
Periodic Fourier interpolation shifts fine samples onto the actual coarse
cell centres before comparing; a separate known-trigonometric-field test
checks this alignment. Successive differences are .0514789 and .0132468,
a ratio of 3.8861. This passes the factor-3 order and .03 finest-difference
gates. These are local planar controls, not convergence of a spherical
pocket or a cosmological calculation.

At N=4096, dt=.0125, halving the axion amplitude gives a scaled waveform
difference of .00016013 after multiplying the smaller field/momentum by
two. It passes the declared .003 linear-amplitude gate. Matched wall recoil
is .00182742 for full amplitude and .000456857 for half amplitude, a ratio
of 3.99997. This is consistent with quadratic weak-field backreaction;
its boundary sensitivity has not yet been established.

## Energy and detector controls

All declared checks pass, including exact initial transfer, exact-zero
axion controls, >98% final coverage, signed exchange and integrated local
continuity, independent output-cadence comparison, order/size and amplitude.

| Fixed-dt spatial case | Returned crossing / incident energy | Final packet coverage |
|---|---:|---:|
| N=2048, dt=.0125 | .99707712 | .99718240 |
| N=4096, dt=.0125 | .99727109 | .99738597 |
| N=8192, dt=.0125 | .99731872 | .99743414 |

The last two returned-crossing measurements differ by 4.76e-5 of incident
energy. The finest window still omits ~.26% of the final packet, which
remains explicitly part of the measurement uncertainty. Integrated axion
continuity residuals are at most 3.43e-5 and output-cadence differences at
most 8.19e-6, within the unchanged 2e-3 gates. Their floor is not evidence
of unresolved time error alone: fixed output cadence and window quadrature
remain in the integrated test.

## Small apparent transmission is unresolved

A read-only endpoint audit is in
`outputs/pockets/planar-small-residual-audit-01/report.json`, produced by
`audit_planar_residual.py`. It partitions the full final axion energy by the
actual dilaton phase, without launching another run.

At fixed dt=.0125 the true-phase energy fractions decrease from 3.511e-6
through 2.120e-7 to 1.313e-8 with successive grid halvings: factors 16.56
and 16.14. At N=4096 their normalized value changes by only .0034% when
the axion amplitude is halved. This is consistent with an amplitude-independent
linear numerical/wall-transient effect whose field amplitude decreases
quadratically with dx. It does not establish the unique cause or exclude
a smaller physical contribution.

These tiny residuals are below the integrated accounting error and cannot
be claimed as a converged physical transmission or long-term leakage rate.
The independent stationary continuum reference predicts full reflection
for the entire initial band. Bulk return and controlled recoil are supported;
a parts-per-million loss measurement is not.

## Remaining gate

The report retains `boundary_control_passed: false` and
`physical_scattering_validated: false`. The second dilaton wall can emit a
small discretisation transient that may reach the returning packet before
the endpoint. Move that wall farther away while preserving target wall,
packet, detector, dx and dt in common physical coordinates; audit any
initial-tail changes. Then compare before interpreting a fully isolated
stationary scattering coefficient or proceeding to moving-wall compression.
Three-dimensional survival, equilibrium scaling, expansion and collision
formation remain unachieved parts of the full goal.
