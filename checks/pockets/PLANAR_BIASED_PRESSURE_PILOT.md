# Field-pressure rebound with positive vacuum bias

2026-10-19. A first native CosmoLattice planar pair shows the coherent axion
field reversing an approaching wall despite continuing inward vacuum
pressure. This is single-grid finite-time evidence; refinement and boundary
controls remain pending. No stable or long-lived pocket is established.

Use the parameters and approximate initial profile specified in
PLANAR_BIASED_PREPARATION.md: canonical fields, Minkowski time, L=512,
N=4096, dx=.125, dt=.0125, T=180, fourth-order spatial Hamiltonian/PV2,
lambda=8, DeltaV=.0002, mFV=.1, mTV=1, fa=10. The identical deterministic
61-mode right-moving packet is used at A=0 and .16. The initial profile is
the boosted degenerate kink, explicitly approximate for positive bias.
The matched empty evolution captures its unforced transient and bias response.

| Quantity | Empty | A=.16 |
|---|---:|---:|
| Final target position | 193.84302 | 237.44581 |
| Final contour velocity | -.388869 | +.160223 |
| Sampled minimum-position time | 180 | 104.5 |
| Duration with velocity >.02 | 0 | 74.5 |
| Position change over final 20 time units | -7.68809 | +3.30988 |

The positive packet satisfies the unchanged reversal criteria (velocity>.02
for >=5 time units and final-20 growth>.5). The empty wall continues inward,
with mean contour acceleration -.0004985, consistent in sign and scale
with the independent initial thin-wall estimate -.0005209. That estimate
is not treated as an exact trajectory reference.

All declared numerical gates pass: exact-zero empty axion, independent
endpoint geometry agreement <1e-12, wall width >=4 cells, matched total
energy closure 5.50e-6 and integrated exchange residual .0001112 relative
to initial packet energy (each below .002). Runtime is 3.985 seconds total
on one thread, below the prospective 15-second cap. No 3D run was launched.

Final axion energy in cells with phi>.9 is 3.34e-6 per transverse area.
This endpoint diagnostic includes penetration and the signed split of
fourth-order gradient energy; it is not cumulative escaped energy and
cannot establish lifetime. A moving-window flux budget will be needed for
a trapped-field slab with repeated interactions.

Evidence: outputs/pockets/planar-biased-pressure-pilot-01 contains the
pre-run contract, native inputs/initial and final fields, full histories,
logs and report with source/build hashes. The shared pressure runner now
accepts --delta-v, defaults to zero, and constructs biased data through the
explicit approximate-profile API. Three biased-seed tests still pass.

Next: separate time/space refinement and displaced-second-wall control at
fixed bias, then distributed trapped-field contraction. A one-sided packet
kick and rebound do not imply an equilibrium pressure-supported pocket.
The full 3D survival, mass-radius scaling and radiation-expansion goals
remain open.
