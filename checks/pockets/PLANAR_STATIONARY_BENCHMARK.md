# Controlled finite-time stationary planar benchmark

2026-10-18. The stationary weak-packet benchmark now passes its declared
implementation, detector, time, space, amplitude, second-wall-sensitivity
and direct-reference gates. This is a finite-time local planar result. It
is not completion of the pocket-survival goal or a precise asymptotic
transmission measurement.

Machine-readable evidence index:
`outputs/pockets/planar-stationary-benchmark-01/report.json`. It names and
hashes the actual reports and their checks, rather than replacing the
original failed pilots or retrospective status flags.

## Second-wall control

`run_planar_boundary_control.py` moves the left antikink from physical x=0
to -256, while the target kink remains at x=256. The domain changes from
[0,512) to [-256,512); dx=.125, dt=.0125, t=0..220, canonical potential,
packet centre/spectrum, and detector remain fixed in physical coordinates.
The compiled backend sees translated coordinates in [0,768), N=6144.
This is not simply rescaling the original geometry.

All four initial field/momentum components match exactly in the common
physical core [32,288]. The original 61 Fourier coefficients are reused;
a quintic taper only suppresses far tails outside that core to make the
larger line periodic. The squared change's free-field energy is 1.102e-17
of initial incident energy. This is an audit of the initial change, not
an assertion that noncompact profiles have strict causal support.

The new packet/empty pair used one thread and 7.34 seconds total, below its
15-second solver budget. Evidence:
`outputs/pockets/planar-second-wall-control-01`, with contract, raw inputs,
outputs and provenance. All declared checks pass. Returned energy differs
from the matched original-domain N=4096 case by -3.86e-9 of initial axion
energy. The final common-core waveform differs by 7.615e-5 in fixed incident
scales, below the declared .003 gate. Measured wall recoil is unchanged.
Final coverage and continuity continue to pass. The exactly zero axion
control remains zero.

This tests one displaced wall and its emitted transient over the chosen
interval. It is not proof of insensitivity to all boundary prescriptions,
late periodic recurrences, or three-dimensional radiation.

## Direct continuum waveform cross-check

`check_planar_reference_convergence.py` uses the saved independent complex
scattering amplitudes and the same 61 incident modes. It reconstructs field
and momentum directly in physical x=[32,224], at least 32 widths left of
the target, where the asymptotic wave expression is controlled. It does
not fit a phase or amplitude to the compiled data. The normalizations are
the fixed initial amplitude .001 and omega0*amplitude=.0005.

At fixed dt=.0125, N=2048/4096/8192 maximum scaled errors are .06887,
.01753,.004290. They improve by factors 3.928 and 4.087 and pass the
predeclared factor-3 and 1% finest-error gates. Evidence and raw comparison
curves: `outputs/pockets/planar-direct-reference-01`. This continuum
reference omits nonlinear recoil; the separate amplitude test bounds the
weak-field discrepancy for these data. It is not a reference for arbitrary
large-amplitude pockets or the right exterior's tiny residual energy.

## What this establishes

The actual coupled field backend reflects a subthreshold weak wave packet,
returns the bulk axion energy, and gives wall displacement with the expected
quadratic amplitude scaling. The finite-time detector accounting and field
waveform have independent, separately refined checks. The returned crossing
on the finest tested grid is .99731872 of incident energy, with .99743414
final coverage; missing detector energy is not declared physical leakage.

The true-phase residual shrinks approximately sixteenfold under grid
halving and lies below integrated accounting error. No parts-per-million
leakage rate, exact reflection coefficient, equilibration or lifetime is
claimed. The next physics benchmark is a moving wall: measure the reflected
frequency shift and laboratory energy exchange, and compare below/above
wall-frame thresholds. Contracting pockets, three-dimensional long-lived
survival, mass-radius scaling, expansion and collision formation still
require their own field experiments.
