# Distributed trapped-field slab: initial-data audit

2026-10-19. Initial data now populate a false-vacuum slab with coherent
standing fields, replacing the one-sided incident packet. No slab evolution
has been run. This is a route toward repeated compression and pressure
support, not evidence of equilibrium or long-lived pockets.

Canonical fields occupy a periodic L=512 line. FV lies between x=128 and
384, with TV outside. Approximate biased kink/antikink profiles have width=1
and inward speed=.05. The potential remains lambda=8, DeltaV=.0002,
mFV=.1,mTV=1,fa=10; time will be Minkowski cosmic time, a=1. Energies are
dimensionless per transverse area. These are coherent classical fields;
there is no isotropy, thermalisation or physical occupation inference.

The axion uses sin(n*pi*(x-128)/256), n=32,36,40, with fixed temporal
phases 0,.7,1.4. Canonical momenta are the corresponding -omega*sin(phase)
quadratures. All modes are odd about the slab centre, so the stress and wall
forcing are symmetric. A smooth quintic taper vanishes within eight wall
widths and reaches unity sixteen units farther inside. These tapered modes
are approximate initial data, not exact finite-wall eigenmodes.

Equal amplitudes are set by sum(A^2*k^2/4)=DeltaV, the untapered linear
phase-average longitudinal pressure estimate. It is not an exact local or
instantaneous pressure-balance condition. In particular the field initially
vanishes near the walls, and pressure reaches them only after propagation.
There is no imposed reflective boundary and no wall force adjustment.

Independent initial-data audits use N=2048,4096,8192 with fixed continuum
coefficients. Energy spatial differences improve by >12 on halving dx.
The finest initial axion energy is .0394486; max|a/fa|=.007073. Energy in
phi>.9 cells is exactly zero initially. Carrier wall-frame penetration
lengths are 1.105..1.176, to be resolved rather than replaced by a Dirichlet
wall. Full Fourier accounting finds ~9.77e-5 of linear FV spectral energy
above the conservative head-on boosted threshold, including taper tails.
This is an initial diagnostic, not a prediction of total escaped energy:
subsequent compression changes frequencies and wall speeds.

Source: planar_distributed_initial.py. Evidence:
outputs/pockets/planar-distributed-seed-02/report.json (adds full tapered
spectrum accounting to the retained seed-01 audit). Three tests check
symmetry, matching empty-wall data, fixed continuum coefficients, pressure
normalization and energy convergence. No compiled evolution is claimed.

Next implement a bounded native pair with the same walls, empty/filled
initial data, and a fixed initial slab-centred energy window. Measure both
wall positions by reflection symmetry, axion exchange and boundary flux;
a passive fixed window cannot by itself measure evolving pocket retention.
Keep prospective energy/exchange/resolution and time/space/isolation gates.
Choose a duration spanning field traversal before claiming repeated wall
interactions; group speeds imply roughly 260 time units per slab crossing.
A short first pilot can only test initial pressure response. Long-time
survival and perturbations, spherical curvature, equilibrium scaling, FRW
and collision formation all remain required by the full goal.
