# Distributed-field slab: first dynamical comparison

2026-10-19. Distributed coherent axion fields reduce contraction of a biased
false-vacuum slab, while both walls continue inward through t=320. This is
initial pressure support, not equilibrium or long-lived pocket survival.

Use PLANAR_DISTRIBUTED_PREPARATION.md for the explicit mode/taper/profile
construction. This pair uses canonical coupled fields, L=512,N=4096,
dx=.125,dt=.0125,T=320 Minkowski time,a=1, order-4 Hamiltonian/PV2,
lambda=8,DeltaV=.0002,mFV=.1,mTV=1,fa=10. Walls start at 128 and 384
with inward speed .05. The filled state contains three coherent tapered
standing modes, with linear phase-average pressure scale DeltaV. The
empty state has identical wall fields/momenta and exactly zero axion.
No particles, thermalisation, isotropy or physical scale inference is used.

| Endpoint quantity | Empty | Filled |
|---|---:|---:|
| Slab length (initial 256) | 163.9224 | 207.9491 |
| Right-wall contour velocity | -.23534 | -.10673 |
| Minimum wall width in cells | 17.09 | 17.50 |

Both endpoint interfaces are measured independently by applying the unique
rising-level diagnostic to the original and reflected fields. Endpoint
reflection/antisymmetry residuals are below 5.6e-13. Historical lengths
2*(right_wall_x-256) are symmetry-based estimates; there is no independent
left-wall history monitor. The minimum sampled lengths equal the endpoints.

Prospective numerical gates pass: matched total-sector closure 5.91e-6,
global axion exchange residual 1.43e-5, and fixed-window axion balance
residual 1.43e-5, each normalized by initial axion energy and below .002.
Independent endpoint geometry, >=4-cell width and zero-axion controls pass.
A stationary passive window centred at 256 has inner radius 144 and outer
160, initially enclosing the field. Net outward axion flux through that
window is 2.64e-5 of initial axion energy, while global axion energy grows
by 16.77% through exchange. The fixed window encloses an increasingly wide
true-vacuum buffer as the walls contract; its small flux is not a moving
pocket leakage fraction or a lifetime measurement.

Runtime is 8.374 seconds total on one thread, below the 15-second cap.
No heavy 3D simulation was run. Evidence:
outputs/pockets/planar-distributed-pilot-01 contains pre-run contract,
input/final fields, native histories, logs and report with source/build
hashes. After both solvers finished, report serialization failed on a NumPy
Boolean. The corrected runner recovered analysis from the existing saved
fields with --analyze-existing; neither solver was rerun. The original
executed runner hash and progress stdout are retained separately.

Next test temporal/spatial refinement and sensitivity to the outer TV
buffer, and add actual moving-boundary retention accounting. T=320 is
about one initial slab traversal, not many repeated interactions. Persistent
contraction is evidence against declaring equilibrium from initial pressure
normalization. Any later turnaround must be followed over multiple crossing
times with leakage/perturbation controls. Spherical survival, equilibrium
mass-radius scaling, radiation expansion and collision formation remain open.

## Coarse time/space sensitivity

Separate coarse comparisons now preserve the same continuum seed, spectrum,
bias, duration and output times. N=2048 at dt=.0125 changes the filled
final length by +.009717 and the maximum right-wall trajectory by .005214.
N=4096 at dt=.025 changes final length by +.0002581 and trajectory by
.0001383. Fixed endpoint field scales [1,.04,1,.02] give maximum differences
.005770 and .001633 respectively. These are diagnostics, not retrospectively
selected pass gates. Both matched empty/filled pairs pass the original
numerical gates and retain slower filled contraction.

The pairs cost 7.84 one-thread solver seconds combined; no fine or 3D
run was launched in this batch. Evidence:
outputs/pockets/planar-distributed-coarse-comparison-01/report.json, with
report hashes for both cases and baseline. compare_distributed_slab.py
checks fixed physics, executable identity and matching history times, and
uses Fourier interpolation to compare cell-centred fields across grids.
Two levels do not establish convergence order. The next controls remain
fine time/space levels and sensitivity to the exterior TV buffer, plus
moving-boundary retention before any survival claim.

A subsequent independent endpoint audit finds ~3.1% outside a window tied
to actual interfaces. See PLANAR_SLAB_RETENTION.md; the small fixed-window
flux above must not be interpreted as this pocket-boundary leakage.
