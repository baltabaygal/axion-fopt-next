# Optional fourth-order planar integration: first waveform-gate pass

2026-10-19. The optional fourth-order spatial Hamiltonian is now implemented
in the actual native NDim=1 CosmoLattice path. The same moving packet passes
the original 1% fixed-speed waveform gate at the same grid and timestep.
This is a new numerical scheme, not a revised physical initial condition or
an overwritten second-order failure. No 3D simulation ran.

## Scoped implementation

`line_laplacian.h` provides a scalar spatial operator for the planar model.
`prepare_line_accuracy_overlay.py` reads the pinned upstream scalar kernel
and writes a user-owned build overlay changing only scalar-Laplacian dispatch.
Potential derivatives remain canonical derivatives. PV2 and other sectors
remain upstream. No vendor file is edited; the generated overlay records
upstream and generated hashes and verifies source content is unchanged.

The existing general CMake line target remains order 2; the dedicated
`build_planar_accuracy.py` creates a separately identified capable build
against the existing CPU libraries. The new build is
`build/pockets-line-accuracy-01/pocket_line_evolve`. Its build.json records
supported orders, compile/link/library/source hashes and pinned vendor
commit. The serial O2 driver build took 6.42 seconds with no fast-math added.
Archived executables are unchanged.

The default `line_spatial_order=2` preserves the previous operator. Order 4
requires the explicit capable build; a normal newly compiled line driver
rejects order 4 without the hook. Runners check build capability metadata
before selecting it. Two ghost cells are allocated for order 4 instead of
one. All initial and evolved state transfers use the framework layout.

`line_spatial_energy.h` and `line_window.h` include the matching signed
hop-2 Hamiltonian and currents. The full-Hessian frequency guard uses
16/(3 dx^2) rather than the order-2 4/dx^2 bound. Default-order arithmetic is
kept on its original branch. Higher-order forces with lower-order diagnostic
energy/current are never presented as a valid configuration.

## Compiled correctness checks

- `planar-order4-backend-validation-01`: actual coupled fields versus independent
  DOP853, N=32,L=8,t=0..0.2, dt=.01/.005/.0025. Finest state error 4.006e-7;
  second-order temporal convergence; energy/source continuity, exact initial
  transfer and periodic translation pass. The changed frequency bound agrees
  with the independent full-Hessian/stencil expression to printed zero error.
- `planar-order4-window-validation-01`: long-bond moving-window diagnostics
  agree with the independent reference within 3.142e-16 of initial window
  energy. Integrated sector residual decreases ~4x per timestep halving to
  1.088e-8. Enabling diagnostics leaves fields and old history values exact.
- `planar-order2-hook-regression-01`: the hook's default order-2 path reproduces
  the old small-run fields and histories exactly; all independent accounting
  and translation gates pass.
- `planar-compiled-spatial-order-01`: six tiny free-mode runs, N=16/32/64,
  L=8, dt=.00025,t=.1,A=1e-6, compare both spatial orders with analytic
  continuum and independent semi-discrete mode solutions. Order-2 error
  reductions are 3.943/3.990; order-4 reductions are 15.674/16.029. The largest
  semi-discrete error is 5.17e-10 of fixed incident scales. Induced dilaton
  displacement is measured at ~1.15e-16, consistent with the weak-field
  approximation. These runs take .061 seconds total on one thread.

These checks exercise the actual operator and ghost layout, not just a
standalone stencil formula. They do not replace physical packet refinements.

## Same-physics moving packet

Evidence: `outputs/pockets/planar-moving-order4-pilot-01/report.json`, raw
inputs, histories, final fields, reference hash and frozen build identity.
N=4096,L=512,dx=.125,dt=.0125,t=0..180,u=.3; unchanged potential, packet
coefficients, amplitude .001, window and independent fixed-speed reference.
A packet/matched-empty pair took 5.03 one-thread solver seconds.

| Observable | New coupled result | Reference |
|---|---:|---:|
| Lab axion energy ratio | 1.84383658 | 1.84392143 |
| Energy-weighted frequency | .93598413 | .93602631 |
| Maximum scaled waveform error | .00108086 | gate .01 |

The direct waveform error decreases from .22506 with order 2 at this same
N/dt to .00108086 with order 4: from 22.5% to .108% in fixed incident
field/momentum scales. No phase, speed or packet amplitude was fitted.
The old order-2 failed gate and all raw outputs remain unchanged.

All declared new pilot gates pass, including the original .01 waveform gate,
1% energy/frequency gates, .002 local continuity/sector closure and >98%
coverage. Axion continuity residual is 5.19e-5 of initial incident energy;
matched dilaton loss versus axion gain closes to 7.25e-6. Final coverage is
99.922%. The empty wall ends at x=202.000888, much closer to the nominal
continuum endpoint 202 than the previous order-2 x=202.093199.

## Remaining work

The pilot retains `moving_scattering_validated: false`. It is a first
same-physics waveform-gate pass, not completion of all fourth-order physical
controls. Time, space, amplitude and second-wall tests must be repeated for
this numerical scheme, with the corrected local final comparison region.
The runners now carry and verify spatial-order identity so a suite cannot
silently mix operators; order-4 spatial improvement requires a factor-12
gate rather than the order-2 factor-3 gate. The fixed-speed waveform threshold
remains .01. Longer contracting pockets, 3D survival, mass-radius scaling,
expansion and collision formation are still unestablished parts of the goal.
