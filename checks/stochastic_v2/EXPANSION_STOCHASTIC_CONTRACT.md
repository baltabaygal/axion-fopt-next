# Expansion followed by stochastic insertion: prospective contract

Declared before implementation/runs, 2026-10-19. Active user-requested goal.
Keep all previous evidence. Use separate owned expanding driver, preserving
the static driver and the polynomial coupled potential unchanged.

## Stage 1: prescribed radiation background

Canonical S,A=f theta, lambda=2, DeltaV=.06, M=.5, f=.2, v=c=1. Potential
U and m_a²=M² S² remain exactly the previous model, not the Marsh exponential.
No temperature dependence, plasma friction, self-consistent Friedmann evolution,
or extra damping. Prescribed radiation domination is a spectator-background
approximation, not a check that the scalar-sector density is subdominant to
an explicitly simulated bath. No Planck/GeV map or reheating claim is made.

Elapsed conformal time u=0 at a=1, initial Hubble H0>0:

    a=1+H0*u; elapsed cosmic time t=u+H0*u²/2
    absolute tau=1/H0+u; absolute cosmic time=1/(2H0)+t
    Hphysical=H0/a²; Hconformal=H0/a
    P_X=a² X'; X'=P_X/a²
    P_X'=a² laplacian(X)-a⁴ V_X, X=S,A.

Snapshots retain the prior binary layout but their time is u and momenta are
P_X, not cosmic velocities: X_dot=P_X/a³. H0 and cosmology are saved alongside
them. Report densities integrated over COMOVING coordinate volume:

    K=integral P²/(2a⁶); G=integral (grad X)²/(2a²); U=integral V
    Ebar=K+G+U (physical total in the box is a³ Ebar)
    Ebar_total'=-Hconformal [6 K_total+2 G_total]
    Ebar_A'=Q_A-Hconformal[6 K_A+2 G_A]
    Q_A=integral (V_A)_S P_S/a².

All rates here are per elapsed conformal time unless explicitly labelled.
Use these balances, never static-energy conservation, for expansion validation.

Prerequisite tests: N4 homogeneous and inhomogeneous coupled states over u0..2
at H0=.2 versus independent DOP853, dt=.02,.01,.005; finest max component
error <2e-5, refinement ratios >3.5, snapshot energy error <1e-11 relative.
Check scale factor/time/momentum transfer, integrated total and sector balances
(finest residual <1e-4 of initial energy), and exact H0=0 regression against
the static backend at the earlier N4 reference endpoints. Save failures.

Single-bubble pilot: R0=8,w=2,theta_i=1, zero initial velocities, periodic
L24,N96,dt=.02,u0..16,H0=.025 (H0/M=.05), a_final=1.4 and elapsed cosmic
t_final=19.2. Controls: N64 at same dt, N96 at dt=.01. Histories .04, slices
.2, initial/final full states, one CPU thread. Energy/source-balance residuals
<.005 of initial energy, independent endpoint energies <1e-11. Temporal final
axion energy change <1%, spatial pilot change <5%; phase fractions reported
without requiring completion. Report comoving versus physical wall scales
and pre-interaction wall resolution; no continuum or settled-vacuum claim.

## Stage 2: stochastic finite-radius seeding

Only proceed after the expansion reference and pilot gates are assessed.
Reuse sample_candidates from src/axion_fopt_next/core.py for the Poisson
candidate process alone. Its Catalogue.wall_speed value is unused metadata;
do not call geometric exclusion or its arrival maps.

Use physical Gamma=H0⁴ exp(beta*t_elapsed) with H0=.025,beta=.4
(beta/H0=16), q(u)=Gamma*a⁴. Positions uniform in comoving volume. Tabulate
q finely and retain the table: sampler exactness applies to its linear
interpolant, not to the original exponential. No extra P_FV factor. Compare
tabulated integrated intensity with independent quadrature (<.001 relative)
and a refined rate table; small seed ensemble checks mean/variance of candidate
counts within declared 5-sigma sampling bounds. Preserve the exact frozen
candidate times/positions for field/timestep/grid comparisons.

At each exact candidate time, evolve the fields to that time, then apply a
declared insertion rule. Require every cell in the seed's physical core
r<=Rseed to satisfy |S|<=.05 (and at least one core cell). This is a deliberate
finite-core, near-FV criterion, stronger than old point exclusion. Record its
measured maximum and rejected events. It is not a thermal bounce acceptance
probability. Multiple candidates at one time are processed in saved order.

For an accepted event, nearest-image physical r=a*d_comoving. Let b be the
tanh seed [1-tanh((r-Rseed)/w)]/2 times a C2 quintic taper from r=Rseed+2w
to r=Rseed+4w, zero outside. Require this support <=a*L/2. Then

    S_new=S_old+(1-S_old)*b; P_S,new=(1-b)*P_S,old
    A and P_A unchanged.

The seed starts at rest relative to the comoving bath (new conformal derivative
zero in its core), not a computed Lorentzian bounce. Its rim may overlap other
tails; no existing converted field is overwritten by a hard wall. Acceptance,
radius, profile, momentum replacement and external energy jumps are part of
the model contract. Log pre/post total/scalar/axion energies for every accepted
event and include their cumulative jumps in the continuity balances. No claim
of self-consistent tunnelling or scalar-sector closed conservation at insertion.

Use exact event alignment via short PV2 steps at event/output times; no rounding
birth times to the field dt and no frozen constant-speed phantom equivalence.
Integration refinement must keep the candidate file fixed. Seed sizes are
physical, so comoving sizes shrink with a. Initial S=0,A=f, both P=0.

Independent insertion tests: small N8,L8,Rseed1.5,w=.5 diagnostic (not physical
supercritical bubbles), fixed events with accepted/rejected/reaccepted cases,
piecewise independent DOP853 with the same insertion map. dt=.02,.01,.005;
finest state error <2e-5, second-order refinement >3, accepted ledger identical,
energy jumps independently verified <1e-10 relative. No-event regression.

Stochastic pilot after these tests: L32,N96,dt=.02,u0..12,H0=.025;
Rseed8,w2, seed 0 (one realization), rate table >=1025 nodes. Controls N64
at same dt and N96,dt=.01, same frozen candidates. Histories .04, slices .2.
Predeclare local trajectory/phase-based acceptance can differ with grid; report
any accepted-catalogue difference as a limitation, not seed equivalence. Energy
balance including expansion and insertion jumps <.005 of max initial/final
energy or summed absolute insertion work; endpoint snapshot accounting 1e-11.
Temporal final E_A difference <2%, spatial <10% if acceptance agrees; otherwise
do not claim a paired physical convergence result. Report nucleation-core and
wall resolution and all parameters. Small pilot only, no abundance, survival,
thermal nucleation or ensemble completion claims. No unbounded scan.
