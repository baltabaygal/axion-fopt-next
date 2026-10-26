# Accessible formation hierarchy: planning audit

2026-10-18. No lattice arrays or simulations were allocated. Evidence:
`outputs/pockets/formation-candidate-audit-04/report.json`. Three small
analytic tests pass in `tests/test_pocket_formation_hierarchy.py`.
Earlier audit versions are retained; audit-04 includes moving penetration,
finite-mode counting, and repeated-encounter estimates.

The coupled backend has passed the small implementation checks in
FORMATION_INTEGRATION.md. This audit asks what physics to test next; it does
not select a demonstrated long-lived pocket or lower any survival gate.

## Candidate family, conditional on later field checks

Canonical excursion units; lambda=8, nominal tanh width w=1, mFV=0.1,
mTV=1, kappa=log(100), fa=10, R0=12, L=48, inward wall speed 0.3, initial
frequency band 0.3..0.5. Compare deltaV=0.01,0.05,0.1 analytically, with
prospective spatial pairs N=192/256, dx=0.25/0.1875. These grids have not
been run. They are prospective finite-thickness field benchmarks, not a
reproduction of the paper's highly relativistic thin-wall hierarchy.

At N=192 the boosted initial tanh 10--90 width is about 8.38 cells.
Wall-frame frequencies at the band edges are 0.4034 and 0.6782, below mTV.
For a locally steady monochromatic evanescent tail, alpha=sqrt(mTV^2-omega'^2)
and x'=gamma(x+ut), so its fixed-lab-time decay length is 1/(gamma alpha).
These lengths span about 4.17..5.19 cells. The shortest initial interior
wavelength spans about 51.3 cells. They pass only initial estimate gates;
actual changing-wall gradients and transmitted fields must be measured.

The four physical field components alone occupy 0.211 GiB at N=192 and
0.500 GiB at N=256. This is a strict storage lower bound, not a process-RSS
estimate or a claim that the full solver fits current free laptop memory.
No production runs should start under the current compute constraint.

The degenerate sigma=1/3 reference gives Rc/R0=8.33,1.67,0.833 across those
biases. Bias/barrier ratios are 0.08,0.4,0.8. Thus thin-wall quantitative
pressure/energy formulae are not controlled for all candidates. Tension,
profile and stress must be measured rather than imported unchanged.

## First-reflection success does not guarantee trapping

At speed 0.3, an incoming lab omega=0.5 has omega'=0.6782 and reflected
lab frequency 0.9219. If it next encounters another approaching wall at
the same speed, its new wall-frame frequency exceeds mTV. The lower band
edge also fails after a small number of such idealised encounters.
The audit records those sequences, stopping at threshold failure.
This is a Lorentz-kinematic estimate for local field scattering, not a
particle simulation or prediction that the coupled wall keeps that speed.
Wall slowing/reversal could change the conclusion and must be observed.

Increasing the exterior mass improves trapping while shortening penetration
lengths and increasing resolution cost. A long run cannot bypass this
tradeoff. First measure reflection/transmission and frequency shifts of
weak field packets against moving dynamical walls, and compare with a
stationary-wall control. Only then test compression and backreaction.

## Finite axion amplitude changes the interior

With deltaV=0.05 and the candidate potential, minimising phi at fixed a
moves the false branch from phi=0 at a=0 to approximately -0.0197 at a=2,
-0.0583 at a=4 and -0.1868 at a=12. The corresponding local mass ratios
relative to unshifted mFV are 0.956,0.874 and 0.650. These are instantaneous
fixed-a extrema, not joint equilibria: dV/da is generally nonzero.

The full Hessian has a negative eigenvalue on the a=12 fixed-a branch,
even though its phi curvature remains positive. The Schur-complement
adiabatic axion curvature agrees with finite differences of the minimised
potential. This is a local curvature diagnostic, not a proven instability
of an oscillating inhomogeneous solution. A seed cannot assume a fixed
interior mass or enforce instantaneous phi relaxation without checking
that approximation's timescale and gradient cost.

A conditional linear, isotropic, monochromatic pressure estimate is saved
in the report to guide amplitudes. It is not an initial-state normalization
or a measured equation of state. For the central bias it gives axion RMS
~1.45 at representative omega=0.4; peaks and mode interference still matter.

## Spectrum and isolation remain gates

A Dirichlet spherical comparison at R=12 contains only ell=1 and ell=2
radial families in the chosen band, with eight angular modes in total.
Finite-wall modes differ. This is not a dense radiation bath: stress
anisotropy must be measured, and complete angular multiplets do not by
themselves guarantee local isotropy or equilibration in one realisation.

The initial interface-to-image gap L-2R0 is 24, whereas the long-lived
threshold t>=20R0 is 240. These are only geometric estimates, not an exact
boundary-arrival guarantee for smooth tails. This box cannot certify
long-term isolation by a simple finite-propagation argument. Boundary
variation or an independently validated outgoing treatment is required;
mass suppression of exterior propagation must not be presumed sufficient.

Next work is field-only initial-profile/scattering preparation at modest
cost. Larger 3D contracting-pocket runs and all survival, scaling, FRW and
collision conclusions remain pending.


The independent continuum field-scattering calculation is now available in
[FORMATION_SCATTERING.md](FORMATION_SCATTERING.md). For the weak-field smooth
degenerate wall it gives ~99.2% transmitted wall-frame flux at the estimated
upper-edge second encounter. This is a reference constraint, not measured
whole-pocket loss or a replacement for dynamical-wall validation.
