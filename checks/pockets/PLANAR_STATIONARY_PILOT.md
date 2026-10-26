# First stationary coupled-wall packet pilot

2026-10-18. One native CosmoLattice NDim=1 weak-packet/matched-empty-wall
pair ran on one thread, within a predeclared 45-second total solver budget.
Actual solver time was 29.20 seconds. No 3D run or parameter scan launched.
Evidence: `outputs/pockets/planar-stationary-pilot-01/report.json`, contract,
input fields/parameters, raw histories, final states, logs and build identity.

Canonical Minkowski fields, N=2048, L=512, dx=.25, dt=.025, t=0..220,
outputs every .5. Degenerate walls lambda=8, deltaV=0, mFV=.1, mTV=1,
kappa=log(100), fa=10. Both fields evolve. The positive-k Gaussian packet
has amplitude 1e-3 and central omega=.5, with the full documented 61-mode
spectrum in PLANAR_SCATTERING_PREPARATION.md. The passive fixed window is
centred at 96, inner radius 80 and outer radius 112. All thresholds below
were declared before the solver runs.

## Observed result and failed gate

Initial axion energy per transverse area is 3.5703723e-6. At the endpoint,
99.99843% of initial axion energy lies in the false-vacuum phase and about
3.512e-6 of initial energy lies in the true-vacuum phase. The latter is an
unconverged small residual, not a validated transmission coefficient.
Total final axion energy is 0.99998781 of initial energy. Almost all energy
has returned to the false-vacuum side, consistent with subthreshold reflection.

The signed current gives outward crossing 0.99997929 and returned crossing
0.93314656, normalized to initial axion energy. Final window energy is
0.93318883. Therefore the declared detector-coverage gate (less than 2% of
axion energy outside the window at the endpoint) FAILS: about 6.68% is
outside. The pilot exits nonzero and its report retains that failed gate.
The threshold is not relaxed and the returned-current fraction is not
reported as a validated reflection coefficient.

Other declared diagnostic gates pass: maximum integrated axion continuity
residual is 4.712e-5 of initial axion energy (<2e-3); independent output-
cadence coarsening changes integrated flux/source by 7.456e-5 (<2e-3).
The integrated window exchange is -1.6665e-6 in the same normalization.
The exact-zero axion control remains zero. These energy gates deliberately
use initial axion energy rather than the much larger total wall energy.

The empty wall remains at x=256. The packet wall ends at x=256.001800,
with positive contour velocity ~2.07e-5, while the empty control's velocity
is numerical zero. This is an observed response with the expected pressure
sign; amplitude/time/space controls are still required to interpret it
quantitatively. Wall-width sampling/relaxation and emitted dilaton waves
remain in both actual field runs.

## Independent finite-band reference comparison

`compare_planar_packet_reference.py` reconstructs the continuum reflected
packet by superposing the actual 61 incident Fourier coefficients and the
complex scattering modes, not just central-frequency energy fractions.
The static, degenerate continuum wall is fixed in this weak-field reference.
All modes lie below the stationary exterior threshold, so its weighted
reflection is one to roundoff. No lattice evolution runs in this comparison.

Evidence: `outputs/pockets/planar-stationary-reference-01`, reference final
field/momentum/gradient/energy and source/pilot hashes. The same endpoint-
shared discrete energy convention is used to compare detector coverage.
The reference itself places only 0.93719479 of initial packet energy in the
chosen window at t=220, though its total packet energy ratio is unity.
Thus detector coverage failure also occurs in an ideal fully reflective
reference; it is not evidence of 6.7% physical loss.

Maximum final field error is 0.06694 of the fixed incident amplitude; momentum
error is 0.06835 of omega0 times that amplitude. This single coarse-grid
comparison does not pass any convergence gate. Long-time lattice dispersion,
wall discretisation/transients and weak recoil have not been separated.
The reference detector fraction differs from the measured one by ~.0040;
neither value licenses altering the original pilot's failed coverage result.

## Next controlled action

Use the finite-band reference to place a detector that actually encloses
the returned packet, keeping the 98% coverage threshold. Then predeclare a
new measurement contract and separately test time, space, amplitude and
boundary effects. Existing raw pilot outputs remain unchanged. A full
scattering validation is still pending. Local planar reflection is only a
prerequisite to contraction/backreaction; it does not establish 3D pocket
survival, mass-radius scaling, radiation-era evolution or collision formation.


A prospective detector redesign now passes the unchanged 98% coverage goal
in the finite-band reference. Offline remeasurement of the saved endpoint
also encloses 99.72%, but there is no corresponding new flux history and
this pilot's failed gate is unchanged. See PLANAR_DETECTOR_DESIGN.md.
