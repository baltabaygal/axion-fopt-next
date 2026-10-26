# Controlled finite-time moving-wall benchmark

2026-10-19. The fourth-order native planar coupled-field benchmark now passes
its declared reference, detector, time, spatial-order, amplitude, energy and
second-wall controls. It tests a weak packet reflecting from an approaching
dynamical dilaton wall. It does not establish a long-lived pocket or finish
the broader field-only programme.

The evidence index is
`outputs/pockets/planar-moving-benchmark-01/report.json`. It lists hashes and
checks of the actual reference, compiled operator/current validations, pilot,
refinement suite and boundary control. Older second-order waveform failures,
initial moving-reference alias failures and broad-region boundary failure
are retained. No thresholds or physical initial data were relaxed.

## Fourth-order physical controls

Raw suite: `outputs/pockets/planar-moving-order4-refinement-01`. Actual
one-thread solver time was 27.07 seconds, within a 35-second wall-time cap.
The physical case remains L=512,t=0..180,u=.3,lambda=8,deltaV=0,mFV=.1,mTV=1,
fa=10, central packet omega=.5,width=16,A=.001. It uses the same fixed-speed
finite-band continuum reference and 1% waveform requirement.

Time at N=4096: dt=.025/.0125/.00625. Successive fixed-scale full-field
differences are .00337094/.000842762, ratio 3.99987. The finest difference
passes .002. This is PV2 temporal convergence; fourth-order *space* does
not imply fourth-order time.

Space at dt=.0125: N=2048/4096/8192. Successive differences are
.0236526/.00152229, ratio 15.5376, passing the declared factor-12 gate.
The comparison preserves the initial continuum mode coefficients and
samples fine fields at coarse cell centres. Direct waveform errors are
.0245893/.00108086/.000861923. The finest error is .0862% of the fixed
incident field/momentum scales and passes the unchanged 1% requirement.
No velocity/phase shift was fitted to the packet.

Half amplitude at N=4096,dt=.0125 gives a scaled axion difference of
.000401662, within .003. All individual pilot accounting, frequency,
energy-gain, detector and zero-control gates pass.

| Fixed-dt spatial case | Lab energy gain | Energy-weighted frequency |
|---|---:|---:|
| N=2048 | 1.84351841 | .93583960 |
| N=4096 | 1.84383658 | .93598413 |
| N=8192 | 1.84385713 | .93599349 |
| Weak-field continuum reference | 1.84392143 | .93602631 |

The reference neglects the finite packet's nonlinear recoil; energy gain
comes from the dilaton. Temporal refinements improve matched-sector energy
closure, and final window coverage remains above 99.9%. The retained
coverage deficit is not silently assigned to irreversible escape.

## Prospectively corrected second-wall control

Raw pair: `outputs/pockets/planar-moving-order4-second-wall-01`. It used
one thread for 6.72 seconds, below its 15-second cap. The left wall moves
from initial physical 0 to -256, while target=256, packet centre=128,
detector centre=88, dx=.125,dt=.0125 and t=180 stay fixed. The larger box
has N=6144,L=768 with physical origin -256.

All four initial components match exactly on [32,288]. The far-tail taper
has squared change energy 1.102e-17 of incident energy. The final comparison
region [62,194] is declared before the run, excluding eight nominal widths
from both known nominal moving-wall paths. It is distinct from initial
transfer overlap and does not include the deliberately swept second wall.

All declared boundary checks pass. Lab energy-gain difference is
1.36e-11; matched recoil difference is -2.84e-14; local full-state
fixed-scale difference is 4.60e-6. Continuity, detector coverage and the
exact-zero control pass. This is a finite-time sensitivity test for one
wall displacement, not a universal boundary or infinite-time isolation proof.
The earlier broad-region failed comparison is unchanged.

## Meaning and next experiment

The verified field equations reproduce the expected moving-wall blueshift
and laboratory energy transfer with backreaction, without particles or a
prescribed wall trajectory in the compiled runs. The actual wall differs
slightly from a strictly constant-speed mirror; weak-amplitude and matched
empty controls bound that difference for this benchmark.

Next increase field pressure in a separately declared finite-amplitude case
and measure wall slowing/reversal against an empty control. That is a
physical hypothesis to test, not an assumed success. Planar turnaround is
not spherical pocket stabilisation. Biased contracting regions, repeated
compression/leakage, long-lived three-dimensional survival, equilibrium
mass-radius scaling, radiation-era expansion and bubble-collision formation
remain separate unfinished requirements of the full goal. No large 3D grid
was allocated during these controls.
