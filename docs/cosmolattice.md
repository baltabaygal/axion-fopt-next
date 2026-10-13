# CosmoLattice adapter feasibility

Target upstream release: v2.0.0, release commit prefix 0dad16e.
Sources checked 2026-10-13:
- https://github.com/cosmolattice/cosmolattice/releases
- https://cosmolattice.net/version13/
- https://github.com/cosmolattice/cosmolattice

The official release confirms new scalar/gauge capabilities and Runge-Kutta
methods. No external stochastic-background API is assumed here. This repository
does not contain a working or compiled CosmoLattice adapter yet.

Prefer transferring one first-arrival field T(x) for common-speed prescribed
bubbles. At every force-evaluation stage, compute the local mass from the actual
stage time and T(x). Updating only at output times or full steps would miss wall
crossings. Higher-order symplectic compositions can query nonmonotonic substeps;
the external field must be queryable without irreversible geometry mutations.

Before writing production code:
1. Inspect v2.0.0 model potential, force and energy evaluators and integrators.
2. Establish the conversion from program time/field/length units to this API.
3. Implement a one-bubble external mass in a reduced scalar model.
4. Check zero mass, homogeneous constant mass, and delayed homogeneous turn-on
   against an independent ODE solution.
5. Refine time and grid spacing across the sharp mass jump; decide whether a
   finite wall profile is needed for reliable convergence.
6. Verify MPI-local indexing, halo layout, axis ordering and GPU memory access.
7. Compare one frozen catalogue with the original lattice simulation before
   enabling stochastic differences.

The NPZ export is a Python interchange artifact, not a format CosmoLattice is
claimed to ingest directly. A C++ reader or agreed conversion format is needed.

## Source inspection and transfer prototype

The pinned v2.0.0 source has now been inspected. The host-side binary reader and
local potential evaluator compile and pass Python/C++ comparisons. The upstream
scalar potential API does not receive stage time, and existing RK cache calls
retain full-step time. See [detailed findings](adapter-prototype.md) before
implementing the actual TempLat-backed adapter. No full solver build or field
evolution benchmark has been completed.
