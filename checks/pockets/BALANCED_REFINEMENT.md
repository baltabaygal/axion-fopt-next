# Balanced field pocket: refinement and box check

Declared 2026-10-17 before launching any of the following runs.
Keep the baseline `balanced-pilot-01` and its failed original gate unchanged.

Question: is the better retention in the balanced-energy pilot numerically
stable, and how much does a larger periodic box change its early evolution?
This is still a static coupled-field study, not a relic lifetime claim.

## Fixed continuum state

Same potential, mode functions, random coefficients and amplitude as
`balanced-pilot-01/contract.json`: R0=20, w=2, lambda=2, DeltaV=0.01,
m_TV=2, 0.15<=omega<=0.6, seed=0, radial h=0.01, r_max=31,
r_cut=30. Positions are relative to each box centre. No gas rescaling
on a new grid or box. Verify radial-mode frequencies and source hashes
against the baseline before evolving. Record all energies and exact transfer.

| Run | N | L | dx | dt | End |
|---|---|---|---|---|---|
| baseline, existing | 128 | 64 | 0.5 | 0.0125 | 80 |
| time | 128 | 64 | 0.5 | 0.00625 | 80 |
| space | 192 | 64 | 1/3 | 0.0125 | 80 |
| box | 192 | 96 | 0.5 | 0.0125 | 80 |

All histories every 0.5, slices every 4, full initial/final states. Same
compiled backend as the baseline, which passes validation-goal-01.
Static time u equals cosmic and conformal time with a=1.

## Gates and interpretation

- Exact initial-state transfer; solver completes; energy drift <5e-4 each.
- Time refinement: maximum absolute radius difference <=0.05 and
  far stored energy fraction difference <=0.005, on shared samples to t=80.
- Spatial consistency: maximum absolute radius difference <=0.25 and
  far stored energy fraction difference <=0.005 over shared samples.
  Two grids test consistency, not continuum convergence.
- Box consistency: maximum absolute radius difference <=0.25 and
  far stored energy fraction difference <=0.005 over shared samples.
  A passing two-box check is not proof of isolated late-time evolution.
- Report original central-row widths unchanged. Also diagnose every .5
  contour edge on the saved plane using local in-plane normals, bracketed
  by nearest .1 and .9 crossings; report missing and ambiguous brackets.
  Missing points cannot be dropped to obtain a pass. The slice omits
  out-of-plane gradients and does not prove 3D resolution. Report a
  diagnostic pass only if all edges are measurable and >=4 cells.

Diagnostic method is in wall_diagnostics.py. Its post-hoc diagnosis of the
baseline finds exterior extra .9 crossings as well as actual thin patches
(minimum 3.30 cells in the plane). Neither the new diagnosis nor a new
refinement can retroactively pass the original pilot.

If spatial or box gates fail, resolve those differences before claiming
long lifetimes. If radius continues to shrink despite early retention,
determine the field energy transfer and wall coupling responsible. No
particle surrogate, imposed reflective boundary, or fixed wall is used.
