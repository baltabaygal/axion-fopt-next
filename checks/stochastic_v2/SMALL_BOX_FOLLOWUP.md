# Small-box follow-up after the L24 pilot

Declared before this follow-up run. The L24,N96,dt=.02 case reached t32 with
F_half=.08067 and no saved full-conversion time. Its energy drift / initial
FV-subtracted excess was .00679, failing .005. The separate dt=.01 control
is retained to evaluate temporal error; do not hide the failed coarse-step gate.

One further smaller-box pilot is authorized by the same request: L=20,N80,
dx=.25, dt=.01, t=0..32, history .04 and slices .2, one CPU thread, 180s
subprocess timeout. Same potential, R0=8,w=2,theta_i=1 and zero velocities.
No damping, altered coupling, clipping or later bubble insertion. Retain the
same global definitions, energy/exchange gates and snapshot-statistics checks
as SMALL_BOX_EXPERIMENT.md. This single follow-up tests whether reducing the
corner distance and changing collision timing suffices to produce any all-cell
true-phase interval; it is not a box-size convergence comparison because the
periodic collision physics deliberately changes with L.

At L20 the closest faces already have S(r=10)=.11920 in the imposed tanh
seed. This appreciable tail and the periodic derivative mismatch are part of
this new initial state and must be reported. It is not an initially isolated
bubble. Do not claim that a further smaller box preserves the old geometry.

If strict all-cell conversion still fails, report it honestly rather than
force S to 1, introduce unrequested damping or conduct an unbounded box scan.
Any need for a dissipative model is a physical decision for discussion. A
transition that leaves field oscillations is different from an irreversible
binary mass-switch model, even if its dominant volume is in the true phase.
