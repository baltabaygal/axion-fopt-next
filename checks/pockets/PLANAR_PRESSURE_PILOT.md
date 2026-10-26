# Finite-amplitude field-pressure pilot: observed wall reversal

2026-10-19. The first finite-amplitude native CosmoLattice planar pilot
observes axion-field pressure reversing an initially approaching dilaton
wall. It is a single-grid observation with pending controls, not a stable
pocket or a long-lived three-dimensional object.

Evidence: `outputs/pockets/planar-pressure-pilot-01/report.json`, predeclared
contract, full input/final fields and trajectory/energy histories. Three
runs use one thread and 6.55 total solver seconds, under a 15-second cap.

## Physical setup

Canonical coupled phi/a fields, NDim=1,L=512,N=4096,dx=.125,dt=.0125,
Minkowski time 0..180, fourth-order spatial Hamiltonian and PV2. Lambda=8,
deltaV=0,mFV=.1,mTV=1,kappa=log(100),fa=10. The target wall starts at x=256
with nominal velocity -.3; the other interface is its periodic partner.
This is a degenerate-vacuum initial-velocity experiment. Contraction here
is not caused by false-vacuum bias or spherical curvature.

The same deterministic 61-mode right-moving packet, central omega=.5,
envelope width 16 and centre x=128, is used at amplitudes 0,.08,.16. There
are no particles, thermalisation assumptions or random phases. Energies
are per transverse area in canonical dimensionless units, with no physical
occupation, mass/radius or cosmological abundance inference. The larger
packet has enough integrated energy to matter to the wall, while its
initial axion angle a/fa remains small; this is finite wall backreaction,
not a claim of strongly anharmonic axion self-interaction.

## Observations

| Incident amplitude | Initial axion energy | Final wall phase velocity | Response |
|---:|---:|---:|---|
| 0 | 0 | -.29914 | Continued contraction |
| .08 | .0228647 | -.13956 | Substantial slowing |
| .16 | .0914586 | +.23521 | Sustained reversal |

For A=.16 the sampled minimum wall position is 227.4272 at t=102.75.
The measured velocity stays above +.02 from t=103.5 through 180, and the
wall advances by 4.7112 during the last 20 units. This passes the declared
reversal criteria: at least five units of positive velocity and more than
.5 units of late position growth. Both position and velocity are used;
a single interpolated sign change is insufficient. The smaller packet
and exact-zero control do not reverse.

The wall is located by its unique rising phi=.5 interface. Canonical
momentum divided by the interpolated gradient supplies a contour-velocity
estimate, with known interpolation bias. Position growth independently
confirms the observed direction. Widths remain above 16.78 cells in all
runs. Ambiguous rising level sets cause an error rather than being omitted.

The fixed-speed weak-field Doppler formula is not used as a success gate
for this nonlinear wall response. Final axion energy ratios are 1.5593 and
1.07039 for A=.08 and .16, different from the weak fixed-speed ~1.844.
The matched empty-wall subtraction shows total-sector closure within
6.34e-6 of initial packet energy; integrated global axion exchange balances
within 1.63e-4, both inside the predeclared .002 gates. The empty axion
field remains exactly zero.

## Diagnostic validation

Optional `line_wall_diagnostics=true` adds wall_x, wall_phase_velocity,
wall_width and wall_width_cells to the line history. It does not alter
forces or prescribe wall motion. The new frozen build is
`build/pockets-line-wall-01`, with source/executable/library hashes and
capability metadata. It compiled serially in 5.34 seconds; old executables
and vendor sources remain unchanged.

The actual C++ geometry agrees with independent Python crossings for static
and boosted walls, including a periodic-cut crossing, and rejects a
multiple-interface fixture. The current diagnostic regression passes with
exact baseline fields/history values when geometry output is off
(`planar-wall-diagnostic-default-regression-01`). With geometry active,
initial/final empty-wall fields and all common long-run energy/history
columns are exactly equal to the old no-trajectory run
(`planar-wall-trajectory-passivity-01`). Endpoint geometry also matches
independent Python results to <=2.85e-14 during the pressure pilot.

## Plot and remaining gates

The inspected exploratory scientific PNG/PDF are in
`outputs/pockets/planar-pressure-plots-02`; plotting code, chart contract
and raw-history hashes are retained. It is a two-panel comparison of
actual position and velocity histories. No fitted trajectory or uncertainty
band is used. Plot-01 is retained; plot-02 moves the annotation away from
a control curve. Matplotlib is used as a scientific fallback because the
optional skill companion/preflight file and Seaborn are unavailable; no
branded-template compliance is claimed.

`turnaround_observed` is true, while `turnaround_validated` and
`three_dimensional_survival_validated` remain false. Next separately refine
time and space, displace the second wall, and check the amplitude response.
Then add a documented vacuum bias and distributed trapped-field initial
state before claiming pressure-supported contracting regions. A directed
packet hitting one planar wall is neither an isotropic pocket gas nor a
mass-radius or lifetime test. The full three-dimensional survival, scaling,
expansion and collision-formation goal remains unfinished.
