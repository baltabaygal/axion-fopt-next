# Optional threaded CPU build for field experiments

Declared before validation or timed comparison, 2026-10-17. The active
stiff-wall suite retains its original Serial executable. Do not replace it
while either paired run is pending. Build the optional OpenMP backend in
a separate directory, using the same pinned upstream sources and strict
floating-point flags. The adapter's default remains Serial, one process,
MPI disabled. No upstream source is changed.

The owned adapter CMake option is AFOPT_KOKKOS_OPENMP, default OFF.
Enabling it explicitly selects Kokkos OpenMP and records the backend in
build-info.json. CPU thread counts are set by OMP_NUM_THREADS, not inferred
from an unverified machine inventory.

## Verification before physical use

1. The complete pocket validation suite must pass with four OpenMP threads,
   including independent ODE comparisons, time refinement, region energy
   and coupled-Hessian guards, and both monitor continuity checks.
2. A shared N=64, L=32, dt=.0125 field state with an inhomogeneous pocket,
   lambda=8, DeltaV=.01, m_TV=2, p=4 must agree between Serial and OpenMP
   to relative 1e-10 in the full final state and 1e-10 in energy histories
   over t=2. Use exact initial-state transfer and hold all other parameters
   fixed. Save source, build, thread count and output provenance.
3. Run the main adapter's existing reference/input checks against the
   optional build too, so the new build option does not silently support
   only the pocket model. Reconfigure the default separately to confirm
   it still selects Serial.
4. Time the same modest evolution at 1, 2 and 4 threads before selecting
   a thread count. Speed is an empirical result, not an assumption. Do not
   compare differently sampled physics runs to claim a speedup.

Failure of these checks keeps the optional backend out of physical use.
The active Serial runs remain the authoritative stiff-wall comparison.
