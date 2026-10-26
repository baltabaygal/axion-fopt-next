# Spectral packet and wall-pair initial data

2026-10-18. `planar_initial.py` supplies deterministic O(N) initial fields;
`prepare_planar_scattering.py` writes inputs and an audit but never starts a
solver. Successful evidence is in
`outputs/pockets/planar-scattering-preparation-02/report.json`, with eight
prepared input directories: stationary/moving walls, N=2048/4096, and
packet/matched empty-wall controls. No new simulation ran. The first
preparation attempt is retained: it wrote inputs but failed while serialising
a NumPy boolean; it is not a successful audit.

## Physical fields and spectrum

The prospective planar domain is periodic, L=512; dx=.25/.125. Degenerate
potential lambda=8, deltaV=0, nominal kink width 1, mFV=.1, mTV=1,
kappa=log(100), fa=10. A kink at x=256 and antikink at x=0 bound the
false-vacuum interval. Periodic image profiles make the endpoints smooth.
For speed u=.3 their widths are Lorentz contracted and the canonical
momenta are -v*dphi/dx, with velocities -u at the right wall and +u at the
left wall. Both momenta are positive as the true phase advances inward.
The stationary case has u=0. These are continuum kink-pair initial data,
not an imposed wall trajectory during the actual coupled evolution.

The packet is centred at x=128. Its analytic Gaussian envelope width is
16, central omega=.5, field amplitude parameter 1e-3. Sixty-one positive-k
Fourier coefficients approximate the envelope with six-sigma truncation.
The real field uses cos(k(x-xc)); its canonical momentum uses
+omega(k) sin(k(x-xc)), with omega(k)=sqrt(k^2+mFV^2). This sign specifies
rightward propagation. The coefficients, physical box and continuum
momenta are held fixed across grids. Retuning momentum to lattice
dispersion would change the initial data and is deliberately avoided.

The full mode band is omega=.1583..8648, including tiny tails. It is broader
than the earlier prospective cavity band; this is a local weak-packet
benchmark, not the candidate cavity gas. There are no random phases or
random seed. It is a classical equation test; quantum occupations and
physical normalization have not been inferred.

## Initial audit

All predeclared preparation checks pass. The continuum left-moving fraction
is below 9e-31. Using lattice dispersion to decompose those same fixed data
gives an unwanted left-moving fraction of 1.008e-7 at dx=.25 and 6.297e-9
at dx=.125. Initial true-vacuum energy contamination is ~2.04e-17; the
fraction within eight nominal widths of either wall is below 2e-18.
These fractions use the separately measured axion energy, avoiding a
subtraction of the much larger wall energy.

The moving-wall packet has a modal initial-energy fraction 5.44e-8 whose
wall-frame frequency is above mTV; the stationary case has zero. This is
an unboosted initial spectral-energy weight, not a wall-frame flux fraction
or measured transmission. The finite bandwidth must be integrated in the
continuum scattering comparison; comparing only the central frequency
would omit it.

The sampled continuum wall is not an exact lattice kink. Maximum
(1-u^2)*laplacian(phi)-U_phi is about .01024/.002641 on the stationary
coarse/fine grids and .01133/.002884 for moving walls. This residual is
recorded rather than counted as a relaxed-wall success. It reduces by
approximately four under refinement. Empty-wall runs must bound resulting
radiation, phase drift and acceleration before attributing differences to
the weak axion field. The tests check boost signs, half-box phase symmetry,
directionality, localization, and common continuum mode coefficients.

## Prospective time windows and remaining gates

Inputs propose dt=.025 and final time 220 for stationary walls, 180 for
moving walls. The packet-centre arrival estimates are 130.64 and 100.02.
A signal from the other wall reaching the target, assuming constant initial
velocities, takes approximately 256 and 196.92. The proposed windows precede
those times. Profiles are noncompact and walls can accelerate, so these
are planning estimates, not a rigorous isolation guarantee. Run a larger
box or matched boundary control before making a quantitative leakage claim.

No prepared input is a completed experiment or a simulation success gate.
Before interpreting reflected/transmitted energy, add local flux or a
validated final-state separation measure, benchmark runtime, evolve the
empty-wall control, refine time separately from space, halve the packet
amplitude to establish the linear limit, and compare the full spectrum with
FORMATION_SCATTERING.md. Conservation refers to total coupled fields;
lab-frame axion gain against a moving wall is not total-energy creation.
The native 1D test can constrain local reflection/backreaction, while 3D
survival, equilibrium scaling, expansion and collision formation remain
separate required stages.
