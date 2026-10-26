# Optional threaded build: numerical agreement and timing

Contract: [THREADED_BACKEND.md](../THREADED_BACKEND.md). The adapter's
default remains Kokkos Serial, one CPU process, MPI disabled. The optional
AFOPT_KOKKOS_OPENMP build uses the same unmodified pinned upstream sources
and strict floating point. It was built separately from the active physics
suite, which was never switched to a different executable mid-run.

Verification completed:

- Accurate Serial and four-thread OpenMP builds each pass 51 pocket checks,
  including independent ODE, Hessian, region-energy and current checks.
- OpenMP passes all 8 existing main-adapter reference checks and its
  invalid-input test. A fresh default configuration and runtime probe select
  Serial with concurrency 1; the optional probe confirms concurrency 1, 2, and 4.
- Shared N=64,L=32,dt=.0125,t=2 nonlinear field runs have bitwise identical
  final states and histories across Serial and OpenMP with 1, 2, and 4 threads.
- Median end-to-end times over three repeats are 0.772, 0.500, and 0.381 seconds
  for OpenMP with 1, 2, and 4 threads. These include matched input/output and were taken
  while another Serial physics run was active. They are not a guaranteed
  large-grid speedup or a production performance claim.
- The compensated-sum probe returns 1 for 1e16+1-1e16. Initial-energy
  replay and static checkpoint continuation checks pass separately.

Raw evidence: `validation-serial-accurate-01`, `validation-omp-02`,
`adapter-omp-validation-01`, `backend-benchmark-01`, `restart-validation-01`
under `outputs/pockets/`. The paired source/executable archives and checksums
are in `threaded-build-01`. No physical survival claim follows from a faster
or numerically matching backend.
