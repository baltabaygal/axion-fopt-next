# Instructions for stochastic nucleation work

Read README.md here, ../AGENTS.md and ../../docs/physics.md before changing
this workflow. This directory owns controlled stochastic-nucleation checks;
it does not replace the homogeneous or single-bubble workspaces.

- Begin with geometry/statistics, using the existing Python sampler. Couple to
  CosmoLattice only after the relevant nucleation checks have passed.
- Keep q=Gamma*a^4 and explicit false-vacuum rejection consistent. Do not insert
  a second P_FV factor into the attempted rate. Distinguish Poisson candidate
  counts from non-Poisson accepted counts.
- Start with a common constant speed and zero initial bubble radius. Changing
  speed histories, bubble-dependent speeds or collision physics requires an
  explicit new physical contract; do not silently reuse phantom equivalence.
- Use an independent analytic/numerical reference. State when the Euclidean
  sphere survival formula is valid and handle finite periodic boxes explicitly.
- Specify units, cosmology, rate, box, time window, speed, seed ensemble and
  percolation/completion definitions in each experiment. Check the candidate
  budget before large runs; expose rate overflow or truncated tails.
- State acceptance criteria before running the comparison. Quantify uncertainty
  across independent realizations and separate probe/grid error from ensemble
  error. Do not treat correlated cells as independent Bernoulli measurements.
- Hold the catalogue fixed for axion lattice/timestep comparisons. For sampler
  or acceleration comparisons, freeze candidate events where possible; a seed
  alone is not a substitute when RNG consumption or rate grids change.
- Preserve the brute-force oracle. Verify accepted events and arrival maps
  before claiming an acceleration is correct; measure both runtime and memory.
- Save parameters, rate tables, RNG/version information, catalogue counts,
  reference calculations, errors and hashes. Use fresh output directories under
  outputs/stochastic/ and preserve earlier evidence.
- Keep initial experiments small and report what each check establishes.
  Large parameter scans, paper abundance/spectrum claims and manuscripts are
  separate tasks, not implicit consequences of opening this workspace.

The user's request establishes this work area. Routine, reversible checks
within a requested experiment do not need an extra approval step. No delegation
is requested by these instructions.
