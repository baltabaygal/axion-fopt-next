# Two-wall retention and finite-time survival screen

2026-10-19. User resumed the field-only programme. Keep one-thread light
compute; no 3D production run. Same continuum distributed seed, box 512,
lambda=8, DeltaV=.0002, mFV=.1, mTV=1, fa=10, initial slab length 256,
inward speed .05. Minkowski time, a=1, energies per transverse area.

Before running: first validate the passive native two-interface monitor
against independent Python endpoint geometry and signed-bond budgets, with
separate dt=.01,.005,.0025 for orders 2 and 4; require agreement <1e-11,
finest integrated residual <1e-5 and refinement ratios >3. Require identical
final fields and old history columns with monitor on/off. Check both sector
continuity equations. Padding 4 uses the existing quintic endpoint window.

Screen at N=2048, dt=.0125, T=1280, cadence=.25 with a 30-second solver cap.
Record independent left/right walls, widths, moving enclosed energies,
exchange, lab flux and boundary advection. Require total drift, global axion
exchange residual and moving axion balance residual <.002 of initial axion
energy; width >=4 cells; endpoint agreement <1e-11. Failures are retained.
No thresholds will be relaxed to obtain a pass. If a wall disappears or
becomes ambiguous the screen fails; saved history must not be called survival.

A rebound is not equilibrium. Signed net moving flux is not irreversible
escape. Outside energy includes penetration and transition-region energy.
T=1280 spans about five initial light-crossing times but is not a 3D lifetime.
A passed exploratory screen still needs fine time/space controls, padding
sensitivity and an exterior-buffer test before physical interpretation.

After the T=1280 baseline, separate N=4096 space, dt=.00625 time, and
L=1024,N=4096 buffer checks all retained the first rebound at sampled
T=716.5 and passed the unchanged gates. Predeclare a T=3840 extension
at N=2048,L=512,dt=.0125 with the same 30-second timeout and gates.
This is a coarse exploratory repeated-interaction test, not a new survival
acceptance gate. Long-duration space/time and exterior controls are required
before accepting any additional oscillations; the earlier buffer test covers
only T=1280. No 3D run is authorized by this experimental contract.
