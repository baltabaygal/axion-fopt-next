# Passive local energy diagnostics in the compiled planar driver

2026-10-18. `line_window.h` is now an optional diagnostic in the actual
CosmoLattice line driver. It does not modify either field's equation or
prescribe wall motion. The separately saved lightweight executable is
`build/pockets-line-window-01/pocket_line_evolve`, with source/build identity
in its build.json. A single changed driver compiled in 2.12 seconds using
existing libraries; previous executables are unchanged.

Parameters: line_window_inner, line_window_outer, line_window_centre,
line_window_velocity. The default outer=0 disables the diagnostic and keeps
the existing CSV columns. Enabled radii obey 0<=inner<outer<L/2. The centre
is c(t)=c0+v*t modulo L. W is one inside inner, zero outside outer and uses
a quintic transition in between. Compact support keeps its periodic
unwrap cut outside the transition. This specified motion belongs to the
measurement window, not to the dynamical dilaton wall.

The ten added history columns are window total/phi/axion energy, corrected
phi/axion flux, signed axion source, laboratory phi/axion flux, and phi/axion
boundary advection. Bond energies are shared symmetrically. Corrected
outward flux is laboratory flux minus advection. Q is positive into the
axion sector: dE_a/dt+F_a=Q; dE_phi/dt+F_phi=-Q. These are net currents,
including incoming fields. They are not automatically irreversible loss
or a reflection/transmission decomposition.

## Small actual-backend validation

`validate_line_window.py` checks current source and executable hashes and
runs the native NDim=1 model with one thread. Raw evidence, endpoint
reference values and provenance:
`outputs/pockets/planar-window-validation-01/report.json`.

N=32, L=8, t=0..0.2, dt=.01/.005/.0025. The compact window has centre=7.7,
velocity=.2, inner=.5, outer=1.5 and crosses the periodic cut in its support.
Initial and final compiled quantities are compared with the independently
implemented `planar_diagnostics.window_budget`. The largest discrepancy is
1.578e-16 of initial window energy. Integrated total/phi/axion continuity
residuals decrease from 1.677e-7 through 4.195e-8 to 1.049e-8, passing the
predeclared 1e-5 finest residual and factor-3 improvement gates. Exchange
and advection are explicitly nonzero, so this is not a vacuous zero-current
test. Endpoint agreement has a predeclared 1e-12 gate.

Full final-field output bytes and every pre-existing history value match
the previously validated no-window runs exactly. Thus the extra diagnostic
is demonstrably passive on these tests. Shifting both fields and window by
seven cells gives identical translated fields and maximum diagnostic change
1.666e-16, within the 1e-12 absolute translation gate. All validation runs
completed in under one second; no full scattering experiment ran.

Next measure a bounded weak-packet case with matched empty-wall control,
using flux signs and arrival windows, finite spectral bandwidth, wall drift,
endpoint residual packet energy and source exchange explicitly. Initial
window energy and elapsed-time slices must be chosen before the physical
comparison. Space/time/amplitude and boundary controls remain necessary.
This supports the field-only pocket programme but does not establish planar
scattering coefficients yet, much less 3D survival or cosmological scaling.
