# Independent planar field-scattering reference

2026-10-18. Computation is a small continuum ODE calculation, not lattice
field evolution. It took under a second and allocated no 3D lattice.
Source: `formation_scattering_reference.py`; reproducible audit:
`audit_formation_scattering.py`. Evidence:
`outputs/pockets/formation-scattering-reference-01/report.json`, with fifteen
complex mode profiles, predeclared tolerances, source hashes and results.

## Background and approximation

Canonical Minkowski fields, lambda=8, deltaV=0, phi=(1+tanh(x))/2, w=1.
At zero axion this is an exact static degenerate dilaton kink. The axion
mass follows the exponential coupling in the formation paper PDF Eq. (2.2):

m^2(x)=0.1^2 exp(log(100) phi(x)), giving mFV=0.1 and mTV=1.

Linearising the axion about a=0 gives, for a=Re[psi(x) exp(-i Omega t)],

psi'' + [Omega^2-m^2(x)] psi=0.

The dilaton response is of order axion amplitude squared and is omitted
in this weak-field reference. This does not prescribe walls in a coupled
pocket run or establish their stability. DeltaV=0 differs from the biased,
accelerating candidate walls; a locally steady boosted approximation needs
its own validation against dynamical field packets.

## Scattering normalization and controls

At the right endpoint impose an outgoing wave when Omega>mTV, or a decaying
evanescent tail when Omega<mTV. Integrate backward with complex DOP853,
then decompose into incident/reflected modes in the left asymptotic region.
Normalize the incident field amplitude to one. Reflection is |r|^2;
transmission is (Re kTV/kFV)|t|^2. These are stationary-wall-frame energy
flux fractions. For a moving wall they are not lab-frame energy fractions:
energy exchange and Doppler factors must be accounted for separately.

Analytic checks pass for a transparent constant-mass medium and a sharp
mass step, below and above threshold. The latter checks complex reflection
amplitude as well as flux. The smooth-wall audit separately varies endpoint
extent 12 to 16 at fixed tolerance and tightens the solver tolerance at fixed
extent 16. All fifteen solutions satisfy the 2e-8 flux balance gate. Complex
reflection changes stay below 1.91e-10 with endpoint variation and 2.92e-11
with tolerance refinement, inside the predeclared 2e-7/2e-8 gates. This is
continuum reference validation, not spatial convergence of a lattice solver.

## Result relevant to repeated compression

| Wall-frame Omega | Reflected flux | Transmitted flux |
|---:|---:|---:|
| 0.403435 | 1.000000 | 0 |
| 0.678208 | 1.000000 | 0 |
| 1.010000 | 0.404993 | 0.595007 |
| 1.254642 | 0.008377 | 0.991623 |
| 1.500000 | 0.000924 | 0.999076 |

For the prospective lab band 0.3..0.5 and approaching-wall speed 0.3,
the first two frequencies are the band edges transformed to the wall frame.
The upper edge reflects to lab omega=0.921911. Encountering another
approaching wall still at speed 0.3 would give wall-frame Omega=1.254642.
Thus the high band edge's second encounter is predominantly transmitted in
this weak-field smooth-wall reference. Merely passing the initial reflection
threshold is not evidence for many-bounce trapping. Near threshold the
reflection fraction remains appreciable, so a binary above/below-threshold
rule alone also loses quantitative information.

This does not predict 99.2% loss of the whole pocket or its spectrum. Actual
walls may slow, reverse or change profile, finite packets have bandwidth and
incidence angles, and finite fields exchange energy with the dilaton. The
next experiment must measure these effects with coupled fields and compare
to this reference in its controlled weak-field limit. Survival, equilibrium,
mass-radius scaling, expansion and formation remain unestablished.
