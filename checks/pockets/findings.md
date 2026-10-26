# Findings: coupled-field relic-pocket programme

Checkpoint: 2026-10-19. Research paused at the user's request. No further
simulations are scheduled by this document. The accepted scientific scope
is [FIELD_GOAL_REVISED.md](FIELD_GOAL_REVISED.md): coupled axion–dilaton
fields only, then controlled survival, equilibrium scaling, radiation-era
expansion and formation from bubble collisions.

## Main conclusion

Coupled classical fields can slow and reverse a contracting planar wall,
including with persistent inward vacuum pressure. Those finite-time
one-sided packet results pass separate time/space refinements and a
second-wall displacement control. Distributed trapped fields also slow
both walls of a contracting slab, but the slab still contracts in the
available evolution. A numerically controlled long-lived three-dimensional
pocket has **not** been established. Neither mass–radius scaling nor
cosmological survival/formation has been demonstrated.

The latest important correction is diagnostic: tiny flux through a large
fixed window does not imply negligible energy outside the pocket. A
wall-following endpoint audit finds about **3.1% of final axion energy
outside the slab window**. This includes penetration and transition-region
energy; irreversible escape and a lifetime cannot be inferred yet.

## Paper interpretation and physical scope

The 2024 formation paper, arXiv:2407.01676, contains a spherical
reflective-particle experiment (Sec. 5, Fig. 4), rather than coupled-field
simulation evidence. The 2026 phenomenology paper does not establish
field-based formation or a damping timescale. Our work tests whether the
mechanism survives field dynamics; it does not reproduce the particle
calculation or infer thermalisation.

The accessible field benchmarks differ substantially from Fig. 4's thin-wall,
highly relativistic hierarchy. Old toy runs used a different mass coupling
and cannot count as reproductions. Physical masses/radii require an explicit
unit map and formation assumptions. No physical pocket mass, radius,
formation efficiency, abundance, or detector rate follows from these runs.
See [model/hierarchy audit](PAPER_FORMATION_AUDIT.md).

## Model and numerical foundation

The current benchmark has canonical fields, with phi=0 in false vacuum
and phi=1 in true vacuum:

    U(phi) = lambda/4 * phi^2(1-phi)^2
             + DeltaV*(1-4phi^3+3phi^4)
    V_ax(phi,a) = fa^2*mFV^2*exp(kappa*phi)*(1-cos(a/fa))
    mTV/mFV = exp(kappa/2).

This quartic family is our explicit benchmark choice, not a uniquely derived
vacuum potential from the paper. Current planar examples use lambda=8,
wall width=1, mFV=.1, mTV=1, kappa=log(100), fa=10. Their coordinates are
Minkowski time with scale factor one; energies are dimensionless per
transverse area. These are coherent fields, not an isotropic thermal gas.

Lattice evidence uses the actual compiled CosmoLattice planar adapter.
Independent Python/ODE references check forces, mixed Hessian, energy
exchange, scattering and spatial/current identities. The optional
fourth-order spatial Hamiltonian has consistent gradient energy, bond
current and frequency bound; PV2 remains second order in time. Vendor
sources are unchanged. Source/build hashes and raw states/histories are
retained in fresh ignored outputs directories.

## Controlled reflection and one-sided pressure

Stationary and approaching-wall weak-field reflection benchmarks pass their
listed finite-time controls. The moving benchmark's finest waveform error
is .0862%, inside its unchanged 1% gate. Its continuum energy-gain prediction
is 1.843921, compared with 1.843857 at N=8192. Temporal/spatial difference
ratios are approximately 4 and 15.54. Earlier failed second-order waveform,
reference-alias and inappropriate boundary-region comparisons remain
archived; they were not converted into passes by relaxing thresholds.
See [moving benchmark](PLANAR_MOVING_BENCHMARK.md) and
[stationary benchmark](PLANAR_STATIONARY_BENCHMARK.md).

One-sided packets use L=512, initial approaching-wall speed=.3, central
omega=.5, envelope width=16, and identical 61 continuum mode coefficients
across grids. Baseline N=4096, dx=.125, dt=.0125, duration=180.

| Vacuum bias | Incident amplitude | Final target velocity | Outcome |
|---|---:|---:|---|
| 0 | 0 | -.29914 | Continues inward |
| 0 | .08 | -.13956 | Slows, no reversal |
| 0 | .16 | +.23521 | Reverses at sampled t=102.75 |
| .0002 | 0 | -.38887 | Accelerates inward |
| .0002 | .16 | +.16022 | Reverses at sampled t=104.5 |

Both strong-packet reversals survive separate temporal/spatial refinements
and displacement of the partner wall. Difference ratios are 4.005/15.326
for zero bias and 4.004/15.074 for positive bias. Turning times stay in the
same sampled bins; finest spatial target trajectories differ from baseline
by at most .000532 and .000602 respectively. Energy/exchange, interface
geometry, wall resolution and exact-zero controls pass. A finite impulse
can reverse a wall without creating an equilibrium pocket.

Evidence: [zero-bias controls](PLANAR_PRESSURE_CONTROLLED.md),
[biased controls](PLANAR_BIASED_PRESSURE_CONTROLLED.md), and raw indexes:

- `outputs/pockets/planar-pressure-controlled-01/report.json`
- `outputs/pockets/planar-biased-pressure-controlled-01/report.json`

## Distributed field: pressure support but continued contraction

A symmetric FV slab initially extends from x=128 to 384 in L=512, with
inward wall speed=.05 and DeltaV=.0002. Three coherent tapered standing
modes n=32,36,40 have fixed phases and canonical momenta. Their untapered
linear phase-average pressure scale equals DeltaV; that is not exact local
pressure balance. The taper initially removes field near the walls, and
these are not exact finite-wall eigenmodes. Its initial spectrum contains
about 9.77e-5 of linear FV energy above the conservative boosted transmission
threshold. Initial axion energy is .0394486 per area.

At N=4096, dx=.125, dt=.0125, t=320:

| Quantity | Empty slab | Filled slab |
|---|---:|---:|
| Initial length | 256 | 256 |
| Final length | 163.9224 | 207.9491 |
| Final right-wall velocity | -.23534 | -.10673 |

Both walls continue inward. This is evidence for pressure support, not
arrest, equilibrium or long-lived survival. The duration covers roughly
one initial slab traversal, not many repeated interactions.

Matched energy closure is 5.91e-6; global exchange and stationary-window
balance residuals are about 1.43e-5 of initial axion energy, below .002.
Endpoint symmetry residuals are below 5.6e-13 and minimum wall widths exceed
17 cells. Both endpoint interfaces are measured independently. Historical
lengths are inferred from symmetry; an independent left-wall history was
not recorded.

The effect persists in separate coarse comparisons. Halving N changes
filled final length by +.009717; doubling dt changes it by +.000258.
These two-level comparisons do **not** establish convergence order; fine
space/time levels and exterior-buffer controls remain outstanding.
The original pair used 8.37 single-thread solver seconds and coarse pairs
7.84 combined. See [distributed pilot](PLANAR_DISTRIBUTED_PILOT.md).

The first report writer failed after both solvers finished because of a
NumPy Boolean serialization error. Analysis was recovered from saved
fields/histories without rerunning the solvers; original executed-source
hash and stdout were retained.

## Retention and leakage: what is actually measured

Global axion energy increases by 16.77% through exchange during the slab
pilot. Net outward flux through the stationary window centred at 256
(inner radius 144, outer 160) is only 2.64e-5 of initial axion energy.
That window includes TV outside the contracting slab and can retain
outgoing energy nearby.

Independent windows following the actual endpoint interfaces instead give:

| Transition padding around interface | Outside fraction of final axion energy |
|---:|---:|
| 2 | 3.056% |
| 4 | 3.118% |
| 8 | 3.543% |

At padding 4, coarse-space/time values are 3.112%/3.120%. The window
includes both translation and changing radius in its time derivative;
finite-difference tests verify that derivative. The outside fraction is
weighted endpoint energy, including penetration and the interface region.
It is not cumulative irreversible escape. The fourth-order local gradient
split is signed, so individual cell contributions also require care.

The missing measurement is a full time-resolved wall-following budget:
energy exchange, lab flux and boundary advection. Two endpoints cannot
separate transient penetration, outward transport and later reabsorption.
See [retention audit](PLANAR_SLAB_RETENTION.md), with evidence at
`outputs/pockets/planar-distributed-endpoint-retention-01/report.json`.

## Remaining work and resume point

1. Complete fine time/space refinement and outer-TV-buffer sensitivity for
   distributed fields, keeping physical initial data and gates fixed.
2. Record time-resolved energy and flux in a window following both actual
   interfaces; independently validate continuity and monitor passivity.
3. Follow compression, any turnaround and leakage over multiple crossing
   times. Test survival rather than accepting a single rebound as success.
4. Establish spherical/3D survival with resolved walls and penetration,
   isolation and perturbation controls, retaining the accepted t>=20 R0
   long-lived gate. Planar evidence cannot establish 3D stability.
5. Only then test dimensionless equilibrium energy/radius relations across
   a controlled family, including stresses, finite tension and isotropy.
6. Extend to radiation-dominated expansion with energy-continuity validation,
   then bubble-collision formation and population evolution. Distinguish
   measured formation outcomes from assumed cosmological parameters.

The user's light-compute constraint remains in force. Recent work uses
short single-thread planar runs and offline diagnostics, with no new heavy
3D production run. Resume from the moving-window retention gap above;
the full scientific goal remains incomplete and is paused.
