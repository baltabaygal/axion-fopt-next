# Memory reduction before further refinement

Declared before implementation, 2026-10-18. Completed N=320/384 runs retain
their original executable and source snapshot. Hardware inventory reports
16 GB RAM. Fine final physical minimum width is .44258, so even a simple
fixed-width estimate needs N>=579 to place four cells across it on L=64;
this estimate is not a prediction of the refined evolution.

Do not launch that grid with the current diagnostic allocations. First
remove redundant cubic storage without changing any fields, forces,
summation order, current definitions or acceptance gates.

First change: replace the tracking monitor's two N^3 window arrays with
four N^2 planes (current/next W and Wdot). Compute each cell's same weights
and preserve face order and compensated reduction order. At N=640 this
reduces this allocation from 4.194 GB to 13.107 MB. Input transfer, full
state extraction and fixed-monitor allocations still need separate audit;
this change alone does not establish capacity for N=640.

Build separately. Before physical use, require all compiled tracking and
structured/default validation gates, plus bitwise field/history agreement
with the frozen previous executable on a moving, nonlinear moderate grid.
Test periodic next-plane wrap, both mass powers and both b values. Record
source/executable hashes and fresh outputs. No new physical convergence,
survival or cosmological claim follows from memory validation.

The plane-buffer implementation passes fresh compiled validation in
`outputs/pockets/validation-lowmem-01`, including independent ODE,
accounting, geometry and second-order continuity checks. A separate N=48
moving pocket crossing the periodic cut has bitwise-identical initial/final
states and histories versus the frozen previous build for all four p=2/4,
b=0/12 combinations (`outputs/pockets/lowmem-agreement-01`). The grid cap
remains 384. Remaining duplicate initial/state and fixed-monitor storage
still prevent claiming safe capacity for N=640 on this host.

Next change: diagnostics read a const host view of the existing four fields
through the same layout mapping instead of constructing a full interleaved
copy. Snapshot serialization uses an 8192-double buffer and slices four
doubles per cell, retaining wire order and headers. Validate all evolved
values before diagnostics just as before. Preserve every loop/reduction
order. Require bitwise initial/final/frame/slice payloads and histories
versus the frozen plane-buffer build, plus the full compiled validation.
Input loading still duplicates the initial state and remains a later step.

This state-view path passes `outputs/pockets/validation-stateview-01`.
`verify_memory_backend.py` records bitwise agreement against the frozen
plane-buffer build for all initial/final/frame/slice binaries and histories
on four moving nonlinear cases (p=2/4,b=0/12) crossing the periodic cut:
`outputs/pockets/stateview-agreement-01`. The modest timing check is slower
with direct indexed reads (roughly twice baseline for most cases; startup
is variable), so no speedup is claimed. At N=640 the removed diagnostic
copy would occupy 8.389 GB. This is an eliminated allocation, not a peak-RSS
measurement or proof that N=640 now fits. The initial transfer vector still
coexists with the model at construction and must be addressed next.

Next input change: map the immutable initial file read-only on POSIX hosts,
retaining exact size, header/grid and finite-value checks before model
initialization. Use memcpy for doubles because the 20-byte header is not
double-aligned. Release the mapping after the existing initialization loop;
keep a buffered fallback elsewhere. Mapping is not a bounded resident-memory
claim: measure peak RSS separately. Require all compiled checks, exact
binary/history agreement and explicit truncated/trailing/nonfinite/header
rejection checks before physical use. The grid cap remains unchanged.

The mapped-input supporting work, started before the revised scientific
goal, finishes with exact binary/history agreement (`mapped-agreement-01`),
compiled validation (`validation-mapped-01`) and six corrupt-input rejection
checks (`mapped-input-rejections-01`) passing. The peak-memory profiling
attempt `mapped-memory-01` does not yield a usable measurement: its baseline
field run completes, but /usr/bin/time -l fails on the sandbox's denied
sysctl kern.clockrate query. No memory-capacity claim or grid increase is
based on this failed profiling attempt. New work follows the accepted
formation-model/hierarchy audit rather than resuming this optimization as
the primary physics objective.
