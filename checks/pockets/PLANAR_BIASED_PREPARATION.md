# Positive-bias planar preparation

2026-10-19. Preparation only: no biased solver evolution has been run.
The next experiment adds a persistent inward vacuum force to the controlled
coherent-packet reversal. It remains a one-wall mechanism test before
constructing a distributed trapped-field slab; it cannot establish pocket
survival or equilibrium.

Initial data use L=512, width=1, speed=.3, packet centre=128, envelope=16,
central omega=.5, A=0 and .16; lambda=8, mFV=.1, mTV=1, fa=10. Introduce
DeltaV=.0002 with the existing formation potential. Fields are canonical,
time is Minkowski cosmic time, a=1, and energies are dimensionless per area.
There are no particles or thermal occupation assumptions.

`planar_biased_initial.py` explicitly samples the boosted degenerate kink
as approximate biased initial data. It retains the strict degenerate-only
guard in planar_initial.py; it does not call this an exact biased traveling
solution. For the target rising FV-to-TV wall,
P=-integral(pi_phi phi_x), so its initial bias contribution is
F=integral(U_bias,phi phi_x)=-DeltaV. Independent quadrature at N=2048,
4096,8192 agrees within 1e-12 relative. Degenerate rest tension is 1/3,
and the thin-wall initial acceleration estimate is -.000520851. This
projection checks sign and normalization, not subsequent dynamics.

Initial false-volume bias energy is .0512 per area; the A=.16 packet has
.0914586 including its consistent fourth-order gradient energy. Evidence:
outputs/pockets/planar-biased-seed-audit-02/report.json. Audit-01 is retained
but its axion-energy entry omitted gradients; use audit-02. Three tests
cover force normalization, explicit approximate-profile construction and
the full initial axion energy.

Before the first solver call, write a fresh run contract with N=4096,
dx=.125, dt=.0125, order=4, T=180, history cadence=.25, one thread and a
15-second aggregate solver cap. Run the empty and A=.16 states with the same
bias and approximate wall profile. Preserve the established .002 matched
energy/exchange gates, unique interface requirement, >=4 width cells and
independent endpoint geometry agreement <1e-12. Observe reversal using the
same velocity>.02 for >=5 time units plus final-20 displacement>.5 criteria;
reversal is a hypothesis outcome, not a numerical-validity gate. Record
leakage and wall acceleration without applying a fixed-speed Doppler fit.
Refinements and displaced-wall controls are required before elevating a
positive pilot outcome. A negative result remains useful evidence and
must not trigger retrospective threshold changes.

Distributed fields, repeated compression, spherical curvature, long-time
survival, dimensionless equilibrium scaling and radiation-dominated
expansion remain outstanding. No 3D run is authorized by this preparation.

## First evolution result

The prospective pair has now run; see PLANAR_BIASED_PRESSURE_PILOT.md.
The preparation-only statement above describes the pre-run state.
