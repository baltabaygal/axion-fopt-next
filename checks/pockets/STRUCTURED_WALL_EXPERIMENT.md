# Separate exterior mass gap from the wall core

Declared before compiled implementation or field runs, 2026-10-17.
The lambda=8 stiff-wall pilot still transfers energy into the dilaton.
Its fixed-sphere export includes translating bound material; a moving
snapshot sphere retains about 97.4% at t=80 and does not yet establish a
radiation rate. Increasing lambda simultaneously narrows
the whole wall and increases its tension and required gas energy. Test a
different, explicit toy dilaton potential that primarily stiffens the
vacuum tails, leaving the central barrier unchanged:

    U_b(phi) = lambda/4 phi^2(1-phi)^2 [1+b(2phi-1)^8]
               + DeltaV (1-4phi^3+3phi^4).

Retain the same canonical real fields and phi^4 axion mass coupling.
This is a new toy potential, not a formula from the reference paper.
Default b=0 must reproduce all existing field equations and trajectories.
The test case is b=12, lambda=2, DeltaV=.01. The power eight preserves
the central barrier and its low-order derivatives, while increasing the
vacuum curvature. For the degenerate potential, b<16 prevents additional
stationary points between the two vacua; the tilted case must also be
checked directly. Do not introduce an intermediate vacuum to obtain a
passing survival result.

Analytically U''(0)=lambda(1+b)/2 and U''(1)=lambda(1+b)/2+12 DeltaV.
An independent quadrature gives sigma=.18730616 and 10-90 thickness
3.9433175 for b=12. The exterior mass is 3.62215 instead of 1.05830
for b=0. Tension rises by 12.4%, rather than the factor two in the
lambda=8 comparison; it is still a recorded confound.

## Preparation and implementation gates

- Verify potential derivatives and Hessian independently, both vacua and
  the number of barriers. Check positivity/boundedness and the default
  b=0 equations. Reject b<0 or b>=16 in this controlled family.
- Construct the degenerate kink by solving in logit(phi), compare its
  tension with independent quadrature and verify spatial residual
  refinement. Do not start with a tanh that has incorrect vacuum tails.
- The compiled implementation must pass the existing validations plus
  new independent ODE, timestep-refinement, energy-region and coupled
  Hessian checks at b=12 before any physics run.
- Preserve old executables and source snapshots; use a separately built
  and validated executable, and record all identities.

## Future field pilot (only after preparation and compiled checks)

R0=20, m_TV=2, phi^4 coupling, a=1, seed=0; Gaussian bound modes with
.1<=omega<=.6, radial h=.01, r_max=31, r_cut=30, solved in the actual
kink profile. Set gas energy to 3(DeltaV+2 sigma/R0)(4 pi/3)R0^3 on
the reference N=192,L=64 grid, then hold amplitude fixed between grids.
Record mode counts and occupations; no thermalisation is imposed.

Use N=192 and 256, L=64, dt=.0125, t=80, monitor radius26,
histories every .25, slices every4 and full states every40. Report both
axion-sector depletion and total net energy leaving the monitor.
The original pilot gates remain: global drift<5e-4, exact transfer,
all wall points measurable and >=4 cells, R>=16 and depletion<=5%.
Monitor continuity/cadence gates are 1e-3. Spatial radius difference
<=.25, far axion fraction difference<=.005 and depletion difference<=.01.
Do not retroactively pass any earlier pilot by changing these criteria.

No short pilot establishes the unchanged ten-crossing-time goal. Isolation,
phase/perturbation dependence, temporal/spatial convergence and a small
total retained-pocket loss rate remain necessary for long survival.
