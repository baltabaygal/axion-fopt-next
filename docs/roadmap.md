# Roadmap and acceptance gates

1. Correctness core (started): sampler, periodic geometry and arrival export.
   Tests provide limited numerical validation; retain the brute-force oracle.
2. Statistics: expanding-RD analytic survival and accepted bubble density;
   rate-grid convergence; percolation threshold explicitly configured; size
   distributions distinguish free-growth radii from collision-limited objects.
3. Acceleration: benchmark spatial indexing against identical candidate events;
   preserve catalogue and arrival results before claiming speed gains.
4. CosmoLattice feasibility: compile actual upstream code and one-bubble adapter;
   pass homogeneous limits and spacetime-resolution convergence.
5. Published benchmark: freeze original commit, cosmology, catalogue, lattice,
   integrator, initial angle, potential and abundance estimator. Compare both
   solvers at multiple resolutions with a predeclared error tolerance.
6. Research: converged density spectra with box-size/UV-cutoff controls,
   ensemble uncertainty, and tests of their survival after the transition.

The supplied Python files are homogeneous ODE/percolation utilities, not the
full published lattice engine. Original lattice source/configs and a frozen
reference output are still needed for gate 5. No paper-reproduction result is
claimed by this scaffold.
