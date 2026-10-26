# Resolved dilaton-mass hierarchy, field-only pilot

Declared before runs, 2026-10-17. Changing the gas band at fixed wall mass
did not produce a candidate. The next test increases the dilaton mass
itself, rather than decreasing the number of occupied gas modes.

Keep the explicit phi^4 model, m_TV=2, R0=20, DeltaV=.01, seed=0,
0.1<=omega<=0.6, Gaussian classical equipartition coefficients. Change
the tanh width to w=1, so lambda=8, m_phi,FV=2 and
m_phi,TV=sqrt(4.12)=2.02978. All initial two-mode source sums are below
the exterior dilaton mass; later nonlinear heating is not excluded.

This also doubles the degenerate-wall tension to sigma=1/3. Record this
confound explicitly. Fix gas energy by the same thin-wall pressure-balance
estimate, Ea=3(DeltaV+2 sigma/R0)(4 pi/3)R0^3, rather than pretending
that equal energy supports both tensions. Neither model is in the bulk
limit R>>R_c: finite tension remains essential to the mass-radius test.

## Paired spatial pilot

Use N=192 and 256, L=64 (dx=1/3 and .25), dt=.0125, t_end=80,
static a=1. The same radial eigenfunctions (h=.01, r_max=31), seed and
amplitude are sampled on both grids. Fix amplitude on N=192, L=64 only;
do not renormalise the finer state. Truncate the negligible axion tail at
r_cut=30. The initial wall 10-90 thickness is 2.197 (6.59 and 8.79 cells).
At omega=.6, exterior axion penetration depth .524 is 1.57 and 2.10 cells.

Histories every .25, slices every 4, full states every 40 plus initial/final.
Enable the independently validated fixed-sphere monitor at radius 26.
All evolution remains the same compiled Hamiltonian field solver. The
monitor is a diagnostic and adds no damping, boundary or wall force.

## Gates and classifications

- Exact initial-state transfer, completion, global energy drift<5e-4.
- Monitor total and sector continuity residuals <=1e-3 of initial monitor
  energy; differences in integrated flux/exchange on twice-coarser saved
  histories <=1e-3. Refine cadence later if measuring smaller losses.
- Report all contour widths and missing brackets in slices and the full
  final 3D state. No unmeasured point is dropped from a diagnostic pass;
  all measured widths must be >=4 cells. Intermediate 3D frames are retained.
- Spatial consistency: maximum radius difference <=.25, maximum far stored
  axion fraction difference <=.005, depletion difference at t=80 <=.01.
  Two grids are a consistency check, not continuum extrapolation.
- A refinement candidate requires R>=16 throughout and axion-sector
  depletion<=5%, plus the numerical and wall diagnostic gates on N=256.
  If inconsistent with N=192, resolve before claiming a controlled lifetime.

No t=80 pilot establishes the t>=400 ten-crossing-time requirement.
Perturbations, phase dependence, time refinement and box/isolation checks
remain mandatory for the goal even if this pilot succeeds. Compare against
the previous lambda=2 runs qualitatively only because tension and bound
spectra change. Report failed gates and continue from the evidence.
