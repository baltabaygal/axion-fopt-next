# Spatial engulfment comparison (before execution)

Keep the pinned all-pairs paper kernel unchanged. Replace only tracer engulfment
with a periodic bubble-to-cell index (CSR), using conservative bounding boxes
and the same float32/nearest-image distance predicate. Fixed 16 cells/axis;
cap at 8 million bubble-cell references, falling back to the original kernel.
No sampling, outer timestep, expansion, rate or RNG changes.

First validate exact surviving indices and order on identical arrays: empty
cases, random configurations, periodic faces/corners, tangency, float32-nextafter
boundary points, large radii and full-box coverage, reordered active indices,
and the capped fallback. Then compare complete catalogues, internal F and probe
arrival maps exactly for all paired runs. No tolerance relaxation for speed.

Paired warm-JIT benchmark: H=beta=1, v=0.5, tau_end=2.5, dt=0.05,
L=3,N=100000; L=6,N=100000 and 800000, seeds 0,1,2. Alternate method order.
Use existing profiling wrappers, exclude projection/checks/output/JIT from
generation timings; include index construction in engulfment. Retain source
snapshots, parameters, hashes, individual timings and catalogues. Each
configuration gets an isolated subprocess and reports process peak RSS including
runtime/warm-up/both implementations; this is not separate per-method memory.
Comparison to earlier timing is secondary; paired measurements are primary.
