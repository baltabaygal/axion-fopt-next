# Prototype: one breathing pocket with axion backreaction

Recorded before any physics run. The question is whether a false-vacuum
region filled with a trapped axion gas survives for many collapse times when
its wall is a dynamical field that feels the axion pressure.

## Model

Two canonical real scalars on a static periodic lattice: the dilaton phi
(vev v=1) and the axion a (unit decay constant). The potential is

```text
V = (1-phi)^2 [lambda/4 phi^2 + dV (1 + 2 phi + 3 phi^2)] + 1/2 m_TV^2 phi^2 a^2
  = lambda/4 phi^2 (1-phi)^2 + dV (1 - 4 phi^3 + 3 phi^4) + 1/2 m_TV^2 phi^2 a^2.
```

False vacuum is phi=0 with V=dV. True vacuum is phi=1 with V=0, so total
energy is the pocket's excess energy. The axion is massless in false vacuum
and has mass m_TV in true vacuum. The backreaction on the wall is the force
m_TV^2 phi a^2 in the dilaton equation; total energy is exactly conserved in
the continuum. The axion is quadratic and free inside the pocket, matching
the free reflective-sphere axion of arXiv:2609.24714, Eq. (3.1).

Scope note (added after the runs started; no parameter, gate or outcome
criterion changed): this prototype starts from a ready-made pocket. It tests
whether axion backreaction can stabilise an existing false-vacuum region,
not whether pockets form from merging bubbles. Reflection of modes with
omega < m_TV makes a bounce likely by construction. The informative results
are therefore the comparisons with the no-gas and leaky controls, leakage
after heating, energy radiated by the wall, asphericity, and which
equilibrium bracket the pocket settles at.

Simplifications: there is no expansion, no nucleation and no thermalising axion
self-interaction. The axion is a classical high-occupation field, not the
paper's thermal particle gas. Pocket radius, wall width, axion wavelength and
m_TV are a reduced hierarchy. Nothing here maps onto physical ARP masses,
radii or transition temperatures.

## Parameters

Units: the false-vacuum dilaton mass sqrt(lambda/2)=1.

| Quantity | Value |
|---|---|
| Lattice | N=128, L=64, dx=0.5, PV2, dt=0.05 |
| Wall | width w=2 (lambda=8/w^2=2), tension sigma=sqrt(lambda/2)/6=1/6, barrier lambda/64 |
| Vacuum energy difference | dV=0.01 |
| Axion true-vacuum mass | m_TV=2 |
| Initial pocket | phi=[1+tanh((r-R0)/w)]/2, R0=20 at the box centre, pi_phi=0 |
| Gas | box Fourier modes with 0.15<=k<=0.6, seed 0, equipartition (field variance 1/k^2, momentum variance 1, independent) |
| Gas window | [1-tanh((r-17)/2)]/2 applied to a and pi_a |
| Gas energy | rho0=0.0464 times (4 pi/3) R0^3, as lattice axion energy including coupling |

rho0 is chosen so that the adiabatic thin-wall equilibrium is R*=17 (0.85 R0).
The thin-wall ODE then predicts a first bounce near R=14.3 at t~18, a
breathing period of about 36 and peak wall speed of about 0.5. The highest gas
frequency then rises from 0.6 to about 0.84, still below m_TV=2.

## Runs

| Suite | Changes from the parameters above | Duration |
|---|---|---|
| control | no gas | t=60 |
| pocket | none | t=300 |
| timestep | dt=0.025 | t=100 |
| spatial | N=96 (dx=2/3), same continuum gas and amplitude | t=100 |
| leaky | m_TV=0.3, gas rescaled to the same lattice axion energy | t=100 |

## Numerical gates (must pass for results to be trusted)

- Validation suite passes (validate.py).
- Initial-state transfer is exact. The driver's dt*frequency guard stays below 0.5.
- Energy drift |E(t)-E(0)|/E(0) < 1e-4 in every run.
- Control: the 3D false-vacuum radius agrees with an independent radial PDE
  reference (dr=0.05, phi=1 held at r=32) to within 0.25 while R_radial >= 8.
- Wall resolution: the 10-90% wall thickness along the central slice row
  stays >= 4 cells.
- Timestep: max |R(dt=0.05)-R(dt=0.025)| < 0.05 for t<=100.
- Spatial: max |R(N=128)-R(N=96)| < 0.5 for t<=100, and their mean radii over
  50<=t<=100 differ by < 0.25. Two grids test consistency, not continuum convergence.

R is the volume-equivalent radius of cells with phi<1/2.

## Physical outcomes (reported either way)

**Relic pocket found** if all three hold in the pocket run:
1. R > 0.4 R0 = 8 throughout 0<=t<=300.
2. The mean R over 200<=t<=300 exceeds 0.5 R0 = 10.
3. Axion energy in deep true vacuum (phi>0.99) is below 1% of the initial
   axion energy at t=300.

**Expected controls:** with no gas, the false vacuum falls below 1% of its
initial volume by t=40. In the leaky run, the axion energy fraction in deep
true vacuum exceeds 50% by t=100, and R(t=100) falls below the pocket run's
R(t=100).

**Recorded without gates:**
- late mean radius compared with the adiabatic (17.0) and energy-conserving brackets;
- breathing period compared with the thin-wall ODE;
- asphericity;
- dilaton fluctuations inside the pocket;
- gas heating (effective wavenumber);
- energy radiated as dilaton waves;
- pocket drift.
