# Globally stiff wall: initial audit only

The declared suite `outputs/pockets/global-stiffness-01` has started its
N=320 field evolution. The N=384 comparison is pending. The executable and
original analysis source remain frozen for this suite.

The independent fresh initial-state audit is saved in
`outputs/pockets/global-stiffness-initial-audit-01`. Exact payload transfer
passes. Slab-wise Python energies match the archived initial calculation
and the compiled energies within 1.60e-16 relative (gate 1e-12). Initial
total energy is 11390.0416918344 and axion energy 7707.37397680696.
The fixed continuum amplitude is 9.19778456652 with 98 angular modes,
seed=0; neither amplitude nor physical parameters may be renormalised on
the finer grid. Full mode labels/frequencies are retained in the suite.

The initial radius is 20.0001684; all 800 initial mid-plane contour edges
have unambiguous 10--90 widths, minimum 5.5224 cells. The initial compiled
dt-frequency bound is .22225, below .5. This is an initial slice check,
not a full three-dimensional resolution or evolution certificate.

The input is a classical Gaussian mode realization, not demonstrated
thermal equilibrium. No completed t=80 run, spatial/time convergence,
perturbation stability, isolation, long lifetime or cosmology result is
proved by these initial gates. Leave the original candidate gates intact.

An independent spectral radial-force audit uses the actual seed/amplitude
and fixed mode actions, without adjusting the ongoing field run. It compares
the Hellmann--Feynman gas derivative against finite radius differences;
relative error is 4.30e-8 (gate 1e-5). Wall inward generalized force is
386.0071 and gas outward force 370.2452, a 4.08% initial inward imbalance.
Minimising only the tanh radius at fixed width/actions gives R=19.7283,
with relative force residual 2.02e-6. This restricted, spherical mean
reference is not a self-consistent wall or a field stability prediction.
The continuum gas energy is 7709.7082, about .0303% above its sampled
lattice value. See `outputs/pockets/global-stiffness-radial-balance-03`.

For comparison, the earlier w=1 case has a 5.19% mean inward imbalance
and restricted equilibrium radius 19.6621
(`outputs/pockets/stiff-wall-radial-balance-01`). Its later field radius
18.78 is not explained solely by this initial mean radial offset. Neither
reference separates non-spherical stresses, energy transfer and radiation.

Early saved mid-plane slices at t=0,4,8,12,16 have unambiguous walls and
minimum widths 5.522,4.640,4.481,4.052,5.033 cells, respectively.
See `outputs/pockets/global-stiffness-early-wall-01`. This sampled early
slice result does not establish intermediate full-3D resolution or pass
the eventual full-duration candidate gate.

## Captured t=20 prefix

`outputs/pockets/global-stiffness-t20-audit-01` saves a complete immutable
history prefix through t=20, its hash, gates and fresh analysis. Early
global energy drift is 1.148e-6; fixed/moving continuity, cadence,
frequency and moving enclosure checks pass. Moving-window stored-energy
retention is .9998565, with minimum inner clearance 1.8023. The radius
decreases to 19.0932; axion-sector energy is .0437% above its initial value.
The t=20 slice is unambiguous, minimum width 4.636 cells. All six saved
early slices pass their diagnostic. This is not settled equilibrium,
full-3D convergence, irreversible loss or long lifetime. The suite keeps
running toward t=80 and then its separately sampled N=384 case.

## Full midpoint snapshot at t=40

The completed full N=320 snapshot is independently audited in
`outputs/pockets/global-stiffness-midpoint-01`. All energy/region components
match the compiled history to 1.56e-16 relative. Global drift through this
time is 2.845e-6. The radius rebounds to 19.5804 after a minimum 19.0716.
Moving-window energy retention is .9992512; axion-sector depletion is
3.2525%, illustrating that exchange between fields is much larger than
the window's net energy change. Minimum false-vacuum inner clearance is
.94854. Neither storage nor net current proves isolated irreversible loss.

All 192736 full-3D contour edges have unambiguous profiles, but 261
(0.1354%) measure below four cells; the minimum is 3.6315 cells, median
5.4604. The declared full-wall resolution check therefore fails at this
coarse-grid midpoint. Do not drop those edges or infer a pass from slice
checks. The pending N=384 case must determine whether finer sampling
resolves the thin patches and reproduces the dynamics. Even if it does,
separate time/box/perturbation and long-duration checks remain necessary.

## Coarse case completes t=80; fine case pending

The frozen suite's N=320 run completes with all nine numerical accounting,
frequency-related evolution guards, current/cadence and enclosure checks
passing. Maximum global energy drift is 1.328e-5. Radius ends at 19.4619,
after a minimum 19.0716. Moving-window retention is .9879558 (1.2044%
net storage reduction); its integrated net outward energy is 136.9728,
with total continuity residual below 1.86e-5 of initial window energy.
Minimum inner clearance remains .59521 throughout saved histories.

The original candidate screen fails: axion-sector depletion is 6.9516%,
above 5%, and wall diagnostics fail. Sector exchange into the axion window
is -532.499 while its net outward axion current integrates to 67.6824;
sector depletion must not be described as whole-pocket loss. The moving
window records essentially no net inward current during this run, but
this alone does not prove absence of periodic-image influence or isolation.
Loss grows late in the history, so do not extrapolate the small initial
loss to t>=400. The N=384 case is now being prepared by the same live suite;
all pairwise convergence and long-lived-pocket requirements remain open.
The final full wall has 186764 measured edges and no ambiguous profiles,
but .7357% are below four cells (minimum 2.496, median 5.485). Several late
saved slices also fail the width threshold. Thus the final thin patches
are more severe than at t=40; their resolution cannot be inferred merely
from the initial width or midpoint result. The fine run must be measured
directly, and may itself require further refinement.

## Late currents and window-support limitations

`outputs/pockets/global-stiffness-current-intervals-01` splits the completed
moving net current into four equal time intervals. Fractions of initial
window energy exported in t=0--20,20--40,40--60,60--80 are .01423%,
.06033%, .20218% and .92583%, respectively. The ending net current is
5.54e-4 of initial window energy per time unit. These are finite-window
net currents, not an irreversible decay law. Their strong late growth
prevents interpreting the modest total loss as a settled small loss rate.
Unresolved wall patches and periodic images remain uncontrolled effects.

A fresh slab-wise support audit in
`outputs/pockets/global-stiffness-window-support-01` reproduces the compiled
weighted energies at t=0,40,80 exactly at reported precision. No cells
with .1<phi<.9 lie outside Rin=26 in these three snapshots. However the
annulus 26<=r<30 contains 6.8303 energy units at t=40 and 66.3470 at t=80;
at t=80 weighted energy exceeds sharp-inner energy by 34.5451 units
(.3033% of initial energy). This quantifies sensitivity to the smooth
window; it does not classify that energy as bound or radiative. Phi<.99
extends far beyond the interface in late snapshots, so a naive .99 contour
cannot define a unique pocket boundary in the presence of waves. The early
enclosure gate tracks false vacuum, not proof of all bound-field support.
