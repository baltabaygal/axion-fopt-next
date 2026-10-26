# Bounded larger-grid support for wall hierarchy checks

The N<=256 prototype cap is a resource guard, not a field equation. A
globally stiffer wall needs smaller dx to distinguish physical retention
from under-resolved wall motion. Extend the controlled cap to N=384,
keeping all time/frequency and input guards. This is not permission to
launch arbitrary large ensembles.

A four-component N=384 double snapshot is 1,811,939,328 payload bytes,
below 2 GiB. Before using this size:

- Stream initial arrays to disk without a second full-size bytes copy.
- Read snapshots directly into their NumPy array, checking header, size,
  finiteness and trailing bytes; do not keep a duplicate raw bytes object.
- Release the Python initial array before starting the compiled process.
- Release the driver's initial transfer vector after the model owns it.
- Verify input transfer by chunked byte comparison, including metadata.
- Keep existing run() behaviour and binary formats unchanged. Check
  roundtrips, corruption detection and equivalence of legacy/prepared runs.
- Keep cases sequential and preserve the current frozen executables. With
  the moving monitor, the compiled working set is expected to be roughly
  4.7 GB at the cap, not including OS/file cache. A memory failure is a
  failed run, not grounds to claim an unresolved lower grid as sufficient.

All new physics runs still require their own predeclared gates. Supporting
N=384 establishes no convergence, stability or cosmological result.

For later full-state wall audits, avoid materialising ndim full gradient
volumes: central differentiation commutes with linear periodic interpolation.
Evaluate the interpolated field one cell either side of each contour point,
then take the central difference. Before using this implementation, compare
sampled gradients with the original full-volume method in 2D/3D (including
points across periodic cuts) to absolute 1e-13 on unit-scale random data,
and require unchanged wall classifications and widths within 1e-11 cells
on representative small profiles. Retain all contour edges, ambiguity
rules, sampling and >=4-cell gates. The current running suite retains its
already-loaded original implementation; any reanalysis must be fresh.

Verification: `outputs/pockets/wall-memory-equivalence-01/report.json`
rechecks the saved N=256 structured-wall final state. All 106810 contour
edges, 4027 ambiguous points, measured widths and the failed classification
are unchanged (reported values identical). Unit tests also compare random
periodic gradients and small curved/missing-threshold profiles. This avoids
ndim full gradient arrays, but does not remove other contour/field buffers
or establish an empirical peak-memory bound for N=384.

`streaming_budget.py` independently sums forward-bond energies and legacy
region splits in x slabs with periodic halos, and hashes the file as a
stream. It retains the old arithmetic FV centroid for exact comparison;
that diagnostic is not a replacement for periodic moving-body retention.
Its temporary arrays scale as slab*N^2 rather than N^3; mapped file pages
and OS cache remain outside this allocation statement. Unit tests cover
both mass powers, structured potential, empty FV, nondividing slabs and
corruption. On the saved N=256 structured-wall final state, all legacy
energy-channel/region values agree to 1.41e-16 relative (gate 1e-12), with
identical file hashes; see `outputs/pockets/streaming-budget-validation-01`.
Future suite analysis uses this implementation and archives it. The
already-running global-stiffness suite retains its loaded original code;
these changes do not alter its frozen contract or evolution.
