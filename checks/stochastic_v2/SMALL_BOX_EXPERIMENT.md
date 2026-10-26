# Periodic single-bubble completion experiment

Declared 2026-10-19 before small-box runs. User requests shrinking the box to
reach full true-phase conversion. Keep the old isolated-bubble evidence intact.
This is now a periodic bubble array: a bubble meets its images. Those collisions
are desired physics here, not an isolated-bubble boundary error.

Keep canonical S and A=f theta, v=c=1, lambda=2, DeltaV=.06, M=.5, f=.2;
same cosine coupling and full reciprocal feedback. R0=8, tanh width w=2,
theta_i=1, zero initial velocities; one imposed seed, no later nucleation.
No expansion, damping, plasma or imposed wall velocity.

Main case L=24, N=96, dx=.25, dt=.02, t=0..32. Coarse spatial control
N=64 (dx=.375), same L,dt,time,seed. Separate time control N96,dt=.01.
History every .04, 2D slices every .2, full 3D initial/final states. One CPU
thread, fresh directories, subprocess guard 180s (240s for the time refinement).
The initial continuum seed has small but non-negligible tails at the closest
faces: S(r=12)=.017986. Keep and quantify them; do not call it an exact bounce
or silently crop/clamp S. The cubic-box corner is r=sqrt(3)*L/2=20.785.
An uncollided sphere estimate alone does not predict completion on this torus.

Add passive global phase diagnostics to the owned compiled driver: minimum,
maximum, mean and variance of S, fractions below S=.5 and below the bare
barrier S=lambda/(2lambda+24DeltaV), and fraction |S-1|<.1. Validate against
independent full snapshots; repeat small-grid DOP853 validation and exact field
regression with the previous compiled driver. No equation or integrator change.

Definitions, declared in advance:

- Operational true phase: S>=.5. F_half=count(S<.5)/N^3. "Full conversion"
  means F_half=0 at a saved global measurement, not a 2D slice or tracer claim.
- Report first full-conversion time, first interval of F_half=0 lasting at
  least 2 time units, endpoint F_half, and any subsequent recrossing. Use the
  actual saved cadence; do not interpolate an exact completion time.
- Also report F_barrier and |S-1|<.1. Being in the true phase is different
  from all points sitting at S=1. Without damping, collisions leave scalar
  oscillations, radiation and a spatially varying axion mass. No frozen
  irreversible conversion or forced constant-mass completion flag is permitted.
- Before collisions, retain the previously defined unique .1/.5/.9 wall
  crossings and >=5-cell-width diagnostic. Once image interaction makes the
  radial-wall description invalid, stop using that width/radius as a single
  spherical-wall observable; global phase fractions remain defined.

Numerical gates: finite fields, exact initial transfer; independent snapshot
energies versus native within 1e-11 of initial total; phase statistics agree
with independently decoded snapshots (fractions exactly, moments within 1e-12).
Energy drift / |initial excess above uniform FV| <.005 and report raw-total
normalization too. Axion exchange residual / initial axion energy <.01.
Time control: final E_A difference <1%, F_half difference <.005; first full
conversion time differs by <=.08 if both have one. Spatial pilot: final E_A
difference <5%, F_half difference <.01, completion time differs by <=.4 if
both have one. These two-grid gates do not prove a spatial continuum limit.
Save failed gates unchanged. If spatial sensitivity fails, use a separately
labelled N128,dx=.1875,dt=.02 run with the same physics rather than relax gates.
If completion has not been reached by t32, record that and prospectively
choose a further bounded extension; do not claim success from geometry alone.
