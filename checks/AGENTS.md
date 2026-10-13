# Instructions for the simple-case checks

This directory is an investigation workspace for controlled scalar-field checks.
Read README.md before changing a case. Current scope includes homogeneous
cosine-potential evolution in radiation domination and the centered periodic
single-bubble setup of arXiv:2402.09501, Sec. IV F (no expansion, quadratic
potential, tanh mass wall). Read single_bubble/README.md for that case.
Stochastic nucleation work lives in stochastic/; read its README.md and
AGENTS.md for the separate statistical validation workflow.

- Preserve the existing solver-independent geometry and adapter contracts.
- Use the actual compiled CosmoLattice adapter for lattice results, and an
  independently implemented analytic/ODE reference for comparisons.
- Distinguish the linear Bessel solution from the nonlinear cosine equation.
  Match both field and derivative at the finite simulation start; theta_i in
  the Bessel formula is the extrapolated regular early-time amplitude.
- State cosmic versus conformal time, scale-factor normalization, mass, box
  length, grid, timestep, initial conditions, and error normalization explicitly.
- Hold physical parameters and timestep fixed when comparing grids. Check
  temporal refinement separately; do not confuse grid independence of a
  homogeneous mode with spatial convergence of an inhomogeneous calculation.
- Check local field/momentum spread, not only spatial means. Normalize errors
  by a fixed amplitude rather than dividing by the oscillating field at zeros.
- Save run parameters, raw histories, provenance, reference curves, tolerances,
  and machine-readable outcomes. Use fresh output directories; do not overwrite
  previous runs or modify third-party source.
- Keep generated solver outputs under the ignored root outputs/ directory.
  Commit-sized summaries and diagnostic figures may live under checks/results/.
- Next cases, only when requested: one bubble nucleating at (L/2,L/2,L/2)
  in a periodic L^3 box, and corresponding non-expanding homogeneous/bubble
  comparisons. Keep stochastic ensembles in their own requested experiments
  under stochastic/. Paper abundance studies remain separate from these checks.

No delegation or extra approval step is required by these instructions.
