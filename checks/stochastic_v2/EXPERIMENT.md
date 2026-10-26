# First seeded dynamical bubble: prospective contract

Declared before the first run, 2026-10-19. Real canonical S and canonical axion
A=f theta, c=1, excursion v=1, Minkowski cosmic time (scale factor=1), periodic
3D box. No stochastic nucleation, thermal bath, friction or expansion.

    V = U(S) + f^2 M^2 S^2 (1-cos(A/f))
    U = lambda/4 S^2(1-S)^2 + DeltaV (1-4S^3+3S^4)
    lambda=2, DeltaV=.06, M=.5, f=.2.

S=0 is a metastable false minimum, S=1 the lower minimum of U. U(0)=.06,
U(1)=0. The bare barrier is at S=lambda/(2 lambda+24 DeltaV)=.367647.
U''(0)=1, U''(1)=1.72. The highest-order coefficient is positive; the
nonnegative axion term does not destabilize large S. For theta=1 the S=0
curvature increases by 2 f^2 M^2(1-cos(1))=.009194; it remains metastable.
The full true minimum at this angle shifts below 1 and must be reported.

The coupled equations are S_tt=lap S-V_S, A_tt=lap A-V_A. No spectator
approximation. The inter-sector source is Q_A=2 f^2 M^2 S S_t(1-cos(A/f));
its integral changes E_A and the opposite source changes E_S. Total energy
includes nearest-neighbor forward-gradient energy compatible with the compiled
CosmoLattice PV2 stencil. No arrival maps or imposed wall velocities.

Initialize one true-vacuum bubble at the centre of L=48, R0=8:
S=[1-tanh((r-R0)/2)]/2, S_t=0, theta=1 everywhere, theta_t=0.
The tanh is a smooth imposed seed, not a computed bounce or exact biased
spherical solution. It has small tails (including at the centre), not compact
support. The reference degenerate planar 10%-90% width is 4.39445; this is
not assumed to remain the evolved width. Degenerate thin-wall tension is
1/6 and estimated critical radius 2 sigma/DeltaV=5.56, only a heuristic for
this biased finite-thickness seed. Growth is measured rather than guaranteed.
The initial seed energy is part of the initial state, not a later source.

Run t=0..6, normal dt=.02, diagnostic/snapshot cadence .5 (choose dt=.025
instead if cadence requires; actual chosen values must be logged). Matched
empty control has theta=0 exactly. Initial angle is fixed across grids.
Spatial pilot grids N=48,64,96; initial reference width spans 4.39,5.86,8.79
cells respectively. These are pilot grids, not a demonstrated continuum limit.
Halve dt separately on N48; add a coarser dt for time-order measurement.
Sources and parameters must be frozen before each run.

Validation gates, declared before looking at outputs:

- Small N4 homogeneous and inhomogeneous states versus independent DOP853,
  dt=.02,.01,.005 over t<=.4: finest max absolute component error <2e-5;
  successive error ratios >3.5; initial transfer exact.
- Independent snapshot energies versus native history: relative error <1e-11.
- Bubble energy drift normalized to initial energy above uniform false-vacuum
  background (absolute magnitude): <.005. Also report raw total-energy drift
  to avoid hiding errors behind the extensive false-vacuum energy.
- Axion-sector integrated exchange residual normalized to its initial energy:
  <.01 with output cadence .1 or finer. Empty axion remains exactly zero.
- Bubble half-height radius from interpolated central axis increases by >.1
  by t=6; this is a demonstration goal, not a numerical validation criterion.
- N48 time-refinement endpoint scaled max-field difference ratio >3.5.
- N64-to-N96 final radius change <5%, axion energy change <10%; report actual
  values and state that two-grid proximity is not spatial-order evidence.
- Finest 10%-90% radial wall spans >=5 cells at all saved times. If an
  oscillatory profile makes the crossing ambiguous, report that rather than
  silently switching contours. Report axis anisotropy.
- Boundary control: compare N64,L48 with N80,L60 (same dx=.75 and same seed)
  at dt=.02; local endpoint field RMS in r<15 changes by <.001 of fixed
  component scales [1,f,1,f], and radius changes by <.5%. This detects both
  image effects and initial-tail sensitivity for this finite time window.

Keep per-process timeouts, one CPU thread, no production runs. A failed gate
remains failed; any revised parameters/criteria require a separately recorded
contract. Save builds, sources, raw states, references, histories and reports
under fresh outputs/stochastic_v2 directories. Compact figures/findings may be
saved here. No relic abundance, terminal speed or continuum claim follows.
