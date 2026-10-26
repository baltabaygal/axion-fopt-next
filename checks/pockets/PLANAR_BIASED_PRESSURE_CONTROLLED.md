# Controlled biased planar rebound

2026-10-19. The positive-bias field-pressure rebound from
PLANAR_BIASED_PRESSURE_PILOT.md persists under independent time and space
refinement and a displaced-second-wall check. This is controlled finite-time
planar rebound; it is not evidence of long-lived pocket equilibrium.

The baseline physical initial data and potential remain fixed: L=512,
T=180 Minkowski time, wall width=1, initial speed=.3, lambda=8,
DeltaV=.0002, mFV=.1, mTV=1, fa=10, coherent packet amplitude=.16 with the
same 61 continuum Fourier coefficients. Every case includes a matched
zero-axion control and uses the compiled order-4 Hamiltonian/PV2 adapter.
Time and space refinements are separate, with fixed comparison scales
[1,.16,1,.08] for phi,a,pi_phi,pi_a and unchanged prospective gates.

| Control | Measured | Requirement |
|---|---:|---:|
| Time differences, dt=.025/.0125/.00625 | ratio 4.00375 | >=3 |
| Space differences, N=2048/4096/8192 | ratio 15.07441 | >=12 |
| Baseline to finest spatial wall-trajectory difference | .000601876 | <.02 |
| Sampled turning-time change | 0 | <.5 |
| Displaced-wall local scaled field difference | 7.90595e-6 | <.003 |
| Displaced-wall target trajectory difference | 1.7053e-13 | <.002 |
| Displaced-wall axion energy ratio difference | 2.2071e-10 | <.0002 |

All runs retain the rebound and pass the original energy/exchange, geometry,
wall-resolution and zero-axion gates. The sampled turning time stays at
104.5. The boundary pair moves the partner wall from physical 0 to -256,
extends the domain to [-256,512], and preserves the initial target/packet
core exactly. Its final comparison interval [80,288] was declared before
running: the biased baseline empty partner reaches about 62.2, so the old
degenerate interval starting at 64 would be inappropriate. This tests one
finite-time boundary displacement, not arbitrary isolation or long-time
recurrence.

Four refinement pairs cost 20.59 one-thread solver seconds under a 30-second
wall cap; the boundary pair costs 6.39 under its 15-second cap. Seven
initial-data/boundary tests pass. No 3D or expansion runs were launched.

Evidence index: outputs/pockets/planar-biased-pressure-controlled-01/report.json.
It hashes the original pilot, refinement and displaced-wall reports, which
retain contracts, raw fields/histories and build/source provenance. Original
pilot statuses are historical and are not overwritten. Shared control
runners now carry the baseline vacuum bias explicitly and verify matching
baseline/refinement provenance before boundary evolution.

Next use distributed trapped fields rather than a single incident packet,
so that repeated compression, pressure support and leakage can be tested.
A finite impulse can reverse a wall even without an equilibrium state;
long-time survival requires a trapped field energy budget and bounded wall
motion over multiple interaction times. Subsequent spherical geometry,
perturbations, mass-radius scaling, FRW and formation remain open.
