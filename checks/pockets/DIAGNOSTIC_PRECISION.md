# Large-grid diagnostic reduction precision

The N=192 stiff-wall run's exact-transfer and ordinary energy-drift gates
pass, but its initial energy-accounting check fails: summing millions of
positive terms sequentially differs from the independent pairwise Python
sum by about 1.5e-12 relative, just beyond the fixed 1e-12 threshold.
This gate remains failed in the original report; do not relax it or
rewrite the original histories to obtain a pass.

Use Neumaier compensated host reductions for subsequent global, regional
and monitor energy/current sums. This changes diagnostics only. The field
potential, force, integrator and timestep guard are unchanged. Never replace
the executable of a pending suite; new binaries are built separately.

Before trusting the correction, retain all usual compiled validation gates,
including independent region energies and sector continuity. Re-evaluate
the exact original N=192 input in a fresh one-step run and compare its t=0
energies with the independent Python calculation at the unchanged 1e-12
tolerance. Check exact input transfer. Record original and corrected errors,
source/build hashes and logs. This one-step diagnostic replay does not
establish survival or retroactively pass the old pilot.

The backend probe also sums 1e16 + 1 - 1e16; the compensated result must
be 1, providing a simple check that finite-only/fast-math optimisation has
not destroyed the compensation. Shared-state Serial/OpenMP tests must
continue to agree after the reduction change.
