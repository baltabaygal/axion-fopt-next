# Completed globally stiff paired-grid test

Both field-only cases in `outputs/pockets/global-stiffness-01` complete
t=80. The raw contract, original source snapshot, initial states and
summaries remain unchanged. Fresh paired analysis is in
`outputs/pockets/global-stiffness-comparison-01/comparison.json`.

| quantity | N=320 | N=384 |
|---|---:|---:|
| final radius | 19.46191 | 19.46892 |
| minimum radius | 19.07161 | 19.06657 |
| moving-window energy retention | 98.79558% | 98.77987% |
| axion-sector depletion | 6.95161% | 6.91407% |
| maximum global energy drift | 1.3283e-5 | 1.2709e-5 |
| minimum final 3D wall width, cells | 2.4960 | 2.6555 |
| final edges below four cells | .7357% | .1624% |

All individual numerical accounting/current/cadence/enclosure gates pass.
All original paired-history gates also pass: maximum radius difference
.009411, maximum far-axion fraction difference .0002464, final depletion
difference .0003754. The additional prerefined-grid moving gates pass:
maximum normalized storage difference .0001840, cumulative-current
difference .0001845, exchange difference .0003591, and periodic centroid
distance .01620. Initial window energies differ by .03454%; the evolved
fields were not renormalised between grids.

Both original candidate screens fail. Sector depletion exceeds 5%, and
both final full-3D wall diagnostics contain widths below four cells despite
having no ambiguous profiles. The fine wall's minimum physical width is
not a simple rescaling of the coarse minimum; further resolution evidence
is required. Passing global histories does not certify these local patches.

This establishes consistent finite-time confinement histories for this toy
field model, not a converged long-lived relic. Sector depletion is not
whole-pocket loss: much of it is transfer to the dilaton. Moving-window net
export grows late, as detailed in the coarse-case audit. Time refinement,
resolved walls, isolation/window support, perturbations and t>=400 survival
remain open, as do mass-radius scaling, expansion and collision formation.
Do not infer physical masses, radii or abundance from this pair.

See [the initial/midpoint/current audits](global-stiffness-initial-01.md)
and [the declared experiment](../GLOBAL_STIFFNESS_EXPERIMENT.md).

## Local contour kinematics

`outputs/pockets/global-stiffness-wall-motion-01` remeasures both final
states with the same contour decisions and records velocities of those
level sets, -pi_phi/|grad phi|. All original widths/classifications are
unchanged. The diagnostic first passes independent translating, contracted
planar-profile tests (speeds 0,.6,.85).

At N=384 the 437 edges below four cells have median absolute level-set
speed .84394, versus .18538 across all measured edges. The minimum-width
point has inferred speed .94297. However 25 thin edges have inferred speeds
>=1, as do some wider edges: these are not accepted material velocities
and are not clipped into Lorentz factors. Motion can correlate with
narrowing, but this does not prove that a boosted kink explains each patch.
Conditional median boost-corrected width is 1.1355 (unboosted planar
reference 1.0986); the unresolved extremes remain uncontrolled.

The fine minimum physical width is .44258. Even holding it fixed would
require N>=579 on L=64 for four cells, and it changes between existing grids.
The host has 16 GB RAM. Further uniform refinement therefore requires
reducing redundant allocations first, rather than merely raising the grid
cap or claiming the current minimum is resolved. See
[the memory contract](../LOW_MEMORY_CONTRACT.md).
