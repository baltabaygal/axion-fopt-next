# Tracer cost baseline (predeclared)

Pinned paper geometry only, H=beta=1, radiation domination, v=0.5,
tau=0..2.5, outer dt=0.05 (50 steps). Boxes L=3 and 6.
Tracer populations 2048,8192,32768,100000, three repetitions with seeds
0,1,2. Same physical analytic F for both volumes: largest radius 1.25 < L/2.
Each configuration runs in a fresh subprocess, loads the pinned declarations,
warms both JIT kernels using 512 tracers, then measures the requested runs.

Instrument class construction, tracer update_state (engulfment), and bubble
evolve_bubbles. Remaining generation wall time includes birth draws, coordinate
gathers, bubble insertion, bookkeeping and catalogue conversion. Independent
arrival measurements at 2048 probe points are timed separately. Exclude
JIT, reference quadrature, file output and plotting from generation timings.
Kernel timers are nested components of generation, not extra additive costs.
Report three individual timings plus medians/ranges, process peak RSS including
Python/Numba and warm-up, and initial resident tracer array bytes (12+8 bytes
per site). Peak RSS is a process high-water mark, not a tracer-only allocation.
Numba engulfment uses 8 threads; bubble radius updates are serial.

Before timed repetitions, require profiling wrappers to reproduce the plain
geometry driver exactly for a warm-up realization (catalogue and internal F).
Save per-realization catalogues, independent probe F, internal F, bubble counts,
phase timings and active-tracer/bubble comparison upper bounds per step.
Require finite nonnegative timers, bounded nonincreasing F, and phase sums
consistent with total generation time. No F accuracy threshold from only three
realizations; this is a cost baseline, not a convergence certification.
Do not change source geometry, RNG ordering, or optimize the reference here.
For subsequent accelerations, freeze the same input events/geometry and compare
to this all-pairs kernel before claiming speedups.

After the first baseline completed, add L=6,N=800000, seeds 0,1,2 with
unchanged gates and timings. This preserves the density of L=3,N=100000
in the eight-times-larger volume. It is an additional cost diagnostic;
existing measurements are retained, not rerun. Run with --density-match-only
into a fresh output directory.
