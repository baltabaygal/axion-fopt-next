# Static field checkpoint continuation

Declared before restart checks, 2026-10-17. For long experiments, reuse
complete saved canonical field/momentum states rather than rerunning an
already completed prefix. This applies only to the static model: its
potential has no explicit time dependence, so a new leg's u=0 can represent
the saved state's original time. It does not justify resetting a background
scale factor or a time-dependent potential.

The PositionVerlet source stores fields and momenta at the same integer
time after every complete drift-kick-drift step. Its sync() only updates a
background derivative; it does not alter static field states. This suggests
that state-only restart is exact, but verify rather than assume it.

Run the same N=6,L=6,dt=.005,lambda=8,DeltaV=.01,m_TV=2,p=4
inhomogeneous field state on the accurate Serial build for t=1. Separately
run to t=.4 on Serial, transfer the complete state, then evolve .6 on the
validated OpenMP build with four threads. Require exact checkpoint transfer
and final-state relative agreement <1e-12. Save hashes, parameters and raw
states. A passing check proves this static checkpoint operation only, not
long-time physical stability or an expanding restart.
