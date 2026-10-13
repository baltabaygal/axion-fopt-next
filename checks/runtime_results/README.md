# 128³ runtime optimization — 2026-10-13

Median executable runtime decreased from **45.423 s** to
**14.165 s**, a **3.207× speedup**
in three alternating before/after pairs. All six final snapshot files match the
original single-bubble reference byte-for-byte. The final particle number is
unchanged; the maximum history-column difference is 2.620e-14, caused
by changing the reduction order for the spatial means.

## Workload and measurement

Centered periodic bubble, L=2*pi, Ngrid=128, v=0.5, quadratic potential,
smooth width=0.15, phi_initial=1, phi_dot_initial=0, no expansion,
double precision. Both executables use dt=pi/320, 1280 steps, t_end=4*pi,
energy/field diagnostics every two steps and initial/final binary snapshots.
The optional mt=2 checkpoint is excluded from both executables.

External wall time includes initialization, evolution, diagnostics and output.
It excludes Python snapshot comparisons and final FFT analysis. Peak RSS is
measured per child process with wait4. The CPU backend is unchanged: Kokkos
Serial on this Apple M5 laptop. No parallel run, reduced output cadence,
float32 conversion, altered timestep, wall thinning or geometry approximation
is used to obtain the speedup.

| Version | Repeat | Wall time, s | Peak executable RSS, MiB |
|---|---|---|---|
| baseline | 1 | 43.753 | 167.11 |
| optimized | 1 | 13.798 | 168.03 |
| optimized | 2 | 14.165 | 168.00 |
| baseline | 2 | 45.730 | 167.12 |
| baseline | 3 | 45.423 | 167.14 |
| optimized | 3 | 14.507 | 168.00 |

The original paper-check run took 44.4 s including Python postprocessing and
an additional checkpoint; this is a different timing definition. Earlier
paired runs also showed host-load variation. The contemporaneous paired
medians above are the reported optimization comparison, not a ratio to a
historical best time.

## Profile breakdown (phase medians, seconds)

| Phase | Baseline | Optimized |
|---|---|---|
| initialization | 0.074 | 0.095 |
| evolution | 25.434 | 11.139 |
| diagnostics | 19.884 | 2.896 |
| snapshots | 0.027 | 0.032 |
| mass update | 17.274 | 0.934 |

The mass-update timing overlaps the other phases; do not add it to the phase
totals. Medians of individual phases need not sum to the median wall time.

## Changes and memory tradeoff

- Evaluate only the selected potential/force/curvature, avoiding sine/cosine
  calls for the harmonic equation.
- Combine moments and energy diagnostics in a single CPU pass with refreshed
  periodic ghosts and the original forward-difference gradient.
- Defer unused endpoint mass updates, while force midpoints, potential queries,
  synchronization and snapshots still obtain the correct mass. Repeated equal
  time requests reuse the current mass field; reversed times are supported.
- Skip tanh calls only where double precision already rounds the result to
  +/-1 (|argument|>=20), retaining the resolved wall and its tails.
- Group exactly equal stored arrival values. This bubble has 14,167 distinct
  arrival times across 2,097,152 cells. Evaluate each mass once per stage and
  scatter it with a uint32 index cache, including padded cells. No radial
  binning, time quantization or rounded arrival key is introduced.

The cache adds roughly 9 MiB of persistent index/level storage at N=128.
Peak RSS is a high-water mark that also includes initialization allocations;
its change need not equal the cache size. `mass_cache=false` disables this
tradeoff, and irregular maps automatically fall back when repetition is
insufficient. This improves runtime at N=128; it does not solve the earlier
1024³ memory constraint.

The first pass without the arrival-level cache reached 1.726× speedup;
[first_pass_summary.json](first_pass_summary.json) retains that measurement.

## Validation and reproduction

The 17 previous checks passed after the optimized cache was added; one new
irregular-map fallback check also passes. The focused bubble suite now has
four checks, including an exact cached/uncached state comparison, independent
DOP853 PDE evolution and energy/work accounting. Full final-build verification
passed all 18 tests with unittest discovery over tests/. The final rebuilt
executable has the same SHA256 as the benchmarked optimized executable.

The paired benchmark compares every final cell (field, momentum, arrival,
mass), all saved history columns and final Fourier particle number with the
original paper-check fixture. [summary.json](summary.json) records timings,
parameters, executable/source hashes, exact snapshot hashes and numerical
comparison errors. Raw six-run outputs are in outputs/runtime-128-cached/.
The original pre-optimization source and binary are archived under
outputs/runtime-128-baseline/. The timing-instrumented baseline changes only
phase timing and timings.json, not the numerical algorithm.

From the repository root, with those archived artifacts and the original
outputs/lee-single-bubble/ reference present:

```sh
PYTHONPATH=src python3 checks/benchmark_runtime.py \
  --baseline outputs/runtime-128-baseline/profiled_evolve \
  --optimized build/cosmolattice/axion_fopt_evolve \
  --output outputs/runtime-128-new
```

Use a fresh output directory. Both solver executables remain single-process
CPU adapters; MPI/GPU acceleration and a float32/two-array redesign have not
been introduced here. Larger-grid or slower-wall speedups require separate
measurements rather than assuming this factor applies universally.
