# Structured vacuum-tail wall: negative survival screen

The fields-only suite declared in `STRUCTURED_WALL_EXPERIMENT.md` completed
on its frozen compiled executable. Raw inputs, source snapshot, build
identity, histories and summaries are under `outputs/pockets/structured-wall-01`.
The potential changes the vacuum tails with b=12 while retaining the central
barrier scale; this is a toy-family test, not the paper's full potential.

| case | final radius at t=80 | axion-sector depletion | relative energy drift |
|---|---:|---:|---:|
| b=0, N=192 | 18.162 | 11.780% | 6.73e-6 |
| b=12, N=192 | 17.429 | 19.392% | 3.24e-5 |
| b=12, N=256 | 17.454 | 19.326% | 3.47e-5 |

All individual numerical transfer, energy and fixed-monitor continuity/cadence
gates pass. Both modified-wall cases fail the original depletion and wall
diagnostic candidate gates. The paired histories pass their declared gates:
maximum radius difference 0.03003, far-axion fraction difference 0.000634,
depletion difference 0.000657. See
`outputs/pockets/structured-wall-comparison-01/comparison.json`.
Agreement of these histories does not rescue the failed wall diagnostic.

A fresh full-3D profile audit (`outputs/pockets/structured-wall-ambiguities-01`)
separates the failure modes without changing thresholds or old reports.
All 120648 initial contour edges are measured unambiguously. At t=40,
1391 lack a .1/.9 bracket within the fixed normal search and 442 have
multiple half-height crossings. At t=80 these counts are 3522 and 505,
respectively, among 106810 contour edges. All measured widths exceed
7.455 cells at t=80. Thus this failed gate is primarily ambiguity of local
profiles, not measured widths below four cells. Missing finite-search
brackets do not prove that a physical wall has vanished; folded profiles
alone do not prove grid under-resolution. They prevent this particular
single-wall measurement from certifying resolution everywhere. No folded
or unmeasured points are silently removed from the failed classification.
Doubling the normal search reach from 8 to 16 units reduces ambiguity from
4027 to 3351 (2809 missing brackets, 542 multiple crossings), while minimum
measured width remains 7.455 cells. The failed classification persists.
This fresh search-sensitivity audit is saved under
`outputs/pockets/structured-wall-search-sensitivity-01`; it does not change
the originally declared reach or acceptance gate.

Axion depletion includes energy exchange with the dilaton and is not total
pocket mass loss. The independent periodic-centroid snapshot audit on N=256
encloses all false-vacuum cells inside radius 26, with final clearance 2.790.
Stored total energies are 4194.8563, 4153.2866 and 3953.8577 at t=0,40,80:
94.255% remains at t=80. See
`outputs/pockets/structured-wall-moving-budget-01/budget.json`.
This is stored energy at three snapshots, not a moving-volume flux,
isolated rest mass, or asymptotic radiation measurement. The fixed sphere
retains 92.799%, illustrating the motion dependence of this measurement.

The tail modification does not pass the intended survival screen. It changes
several physical properties together and cannot establish a unique leakage
mechanism. The next experiment must be separately declared and should
address the nonlinear wall response and resolution rather than infer a
cosmological relic from this negative result. No t>=400 lifetime,
perturbation stability, mass-radius scaling, expansion or formation claim
is established.

Supporting numerical work: the moving weighted-volume monitor passes fresh
compiled validation in `outputs/pockets/validation-tracking-02` after the
streamed I/O/resource changes. Full repository tests with `PYTHONPATH=src`
pass (71 tests, 11 optional skips). Larger-grid capacity has no physical
convergence claim attached.
