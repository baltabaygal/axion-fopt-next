# Field-only early evolution: refinement and preparation

Contracts: [BALANCED_REFINEMENT.md](../BALANCED_REFINEMENT.md) and
[RELAXED_PILOT.md](../RELAXED_PILOT.md), recorded before their field runs.
All runs use the compiled two-real-field model, unchanged phi^2 coupling,
static a=1, lambda=2, DeltaV=.01, m_TV=2, initial radius about 20.
Raw outputs and source/build identities are in `outputs/pockets/` under
the directories named below. No particle implementation is used.

## Fixed continuum-state comparisons to t=80

Baseline: `balanced-pilot-01`, N=128, L=64, dt=.0125.
Refinement suite: `balanced-refinement-01`; every run holds mode functions,
seed and amplitude fixed. No grid-dependent gas energy normalisation.

| Change | Maximum radius difference | Maximum far stored axion energy difference / initial E_a |
|---|---:|---:|
| dt=.00625 | 0.000244 | 0.00000515 |
| N=192, L=64 | 0.06152 | 0.003147 |
| N=192, L=96 (same dx=.5) | 0.007730 | 0.000567 |

All pass the predeclared early history consistency thresholds; total energy
drift is <=1.02e-5. This is two-grid/two-box consistency, not continuum or
isolated lifetime proof. The new wall diagnostic is a separate check:
all suites have some ambiguous slice contour brackets. The original failed
central-row gate is not retroactively changed.

## Wall diagnosis

The central row has two .5 crossings throughout the balanced baseline;
extra .9 crossings arise in exterior ripples. The post-hoc all-contour
normal check (`wall-audit-01`) also finds actual thin patches: minimum 3.30
cells in the plane. The final full-3D check (`wall-audit-3d-01`) measures
23244 contour-edge points; 352 lack an unambiguous .1-.9 bracket and the
minimum among the measured points is 3.04 cells. Missing points are not
discarded from pass decisions. Thus, the central slice alone misses worse
wall-resolution regions. No wall-resolution pass is claimed.

## Energy channels

For the balanced baseline, axion-sector energy (including the coupling
potential) falls from 2680.83 to 2027.90 at t=80. Dilaton-sector energy
increases by the same amount within the conserved total. Only 85.15 of
axion energy is in the geometric far field at the final time. Far dilaton
energy grows from 0.029 to 693.00, using r>R_eff+6. The larger L=96 box has
738.91 of far dilaton energy. See `field-budget-01`.

These are stored-energy allocations at two times, not cumulative fluxes.
They show why the far axion fraction alone misses the dominant energy
transfer. The field coupling allows axion energy to enter dilaton waves.
The detailed frequency-transfer mechanism has not yet been isolated.

## Mean-wall preparation does not fix the loss

An independent radial fixed-action ensemble reference converges at h=.1
and .05 to R_half=19.65857 and 19.65880. Energies differ by 0.0041% and
force-density residuals after polishing are <8e-14. Directional energy/force
checks pass. This is an equilibrium reference, not evolved field evidence.

The `relaxed-pilot-01` compiled pair uses the same per-harmonic actions and
random phases, with tanh versus self-consistent mean-wall profiles. It is a
new fixed-action classical ensemble, not the previous Gaussian realisation.

| Initial wall | Initial radius | Final radius at t=80 | Axion-sector depletion | Final far axion storage / E_a(0) |
|---|---:|---:|---:|---:|
| tanh | 20.000 | 16.666 | 23.44% | 3.37% |
| relaxed mean wall | 19.670 | 16.602 | 24.79% | 3.09% |

Both conserve total energy to <1e-5 and complete, but both fail the slice
wall diagnostic. The predeclared improvement criterion fails. Mean-wall
relaxation alone does not prevent the transfer and shrinkage in this model.

## Implication

Do not spend a long lifetime run on this state or infer physical mass/radius
from it. Next, compare a documented mass-function family with vanishing
false-vacuum curvature (phi^4 coupling) against the phi^2 baseline, validating
the changed force and compiled energy accounting independently first.
This tests a structural source of coupling rather than tuning a lifetime.
No long-lived relic, cosmological abundance, or expanding survival result
has yet been established.
