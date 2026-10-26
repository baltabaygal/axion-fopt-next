# Axion/dilaton frequency hierarchy, field-only pilot

Declared before runs, 2026-10-17. The broad-band phi^4 pilot still loses
13% of axion energy to the dilaton by t=80. Test whether the available
quadratic axion source frequencies influence this transfer.

The exterior dilaton mass is sqrt(lambda/2+12 DeltaV)=1.05830, obtained
by differentiating the explicit model potential at phi=1. The phi force
contains a^2, so two axion modes can drive sums/differences of their
frequencies. This motivates a threshold test, not an assumption that the
initial band stays fixed during nonlinear evolution. No claim that this
channel is the only source of radiation is made.

## Fixed parameters and deliberately changed state

Use the validated compiled phi^4 model: lambda=2, DeltaV=.01, m_TV=2,
R0=20, tanh width w=2, wall initially at rest. N=128, L=64, dx=.5,
dt=.0125, a=1, t_end=80 (two initial diameter crossings).

Bound-mode Gaussian field and momentum coefficients use seed=0,
equipartition variance 1/omega^2 and 1 respectively, high-occupation
classical interpretation only. Solve radial functions with h=.01,
r_max=31; truncate the exponentially small field tail at r_cut=30.
The initial gas energy is 2680.825731 on the reference N=128, L=64
grid. Its amplitude is held fixed if the same state is later sampled
on another grid or in another box.

| Case | Initial omega band | Largest initial pair sum |
|---|---|---|
| soft | .1 to .35 | .70, below exterior dilaton mass |
| intermediate | .1 to .45 | .90, below exterior dilaton mass |
| broad | .1 to .60 | 1.20, above exterior dilaton mass |

The lower cutoff is .1 to retain the fundamental spherical mode. In the
older phi^4 pilot, the .15 cutoff excludes the fundamental near .1497;
therefore this broad control is a new initial state, not a rerun of the
old pilot. State counts, actual eigenfrequencies and amplitudes must be
recorded. Different bands and normalisation change the number of modes,
amplitude and angular pressure fluctuations; this is an explicit confound.
Do not infer a single-channel cause from a band comparison alone.

## Gates and report

Exact initial-state transfer; completion; energy drift<5e-4.
Record slice contour and final full-3D wall diagnostics, counting every
missing/ambiguous bracket as a diagnostic failure rather than omitting it.
All measured 10-90 wall thicknesses must be >=4 cells for a diagnostic
pass. Slice/final samples do not prove resolution at unsaved times.

Compare radius loss, total axion-sector depletion, stored far axion energy,
far dilaton storage, asphericity and characteristic gas frequency. A soft
case is a candidate for refinement only if R>=16 throughout, depletion at
t=80 <=5%, and it passes the numerical and wall diagnostics. Report partial
improvement and failed gates either way. No t=80 pilot establishes the
ten-crossing-time survival requirement. Refine any useful candidate,
perturb its state, test box dependence, and extend duration separately.
