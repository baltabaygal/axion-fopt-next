# Additional moving-current grid comparison

Declared before the global-stiffness N=384 state/evolution is available,
2026-10-17. The N=320 case is already complete and fails its original
candidate screen. These additive gates do not change that result or replace
any original radius, depletion, wall or numerical gate.

At shared saved times, compare each run's moving weighted energy divided
by its own initial window energy. Maximum difference must be <=1e-3.
Compare cumulative trapezoidal total net outward current and cumulative
coupling exchange, likewise normalized to each initial window energy;
each maximum difference must be <=1e-3. For L=64, maximum periodic
distance between phase-weighted centroids must be <=.25. Record raw initial
energy differences as well; normalization is not energy renormalisation
of the evolved fields. The suite already uses the same continuum amplitude.

Only compare energies, depletion and currents at the same times. A prefix
may be reported explicitly but fails the full-duration gate. Two grids
test consistency; they do not supply a continuum extrapolation, repair
under-resolved walls, establish temporal convergence or certify isolation.
These 1e-3 absolute energy-fraction thresholds are numerical error targets,
not an allowed physical decay rate. An eventual small-loss claim must
resolve its measured rate relative to discretisation/current errors.
