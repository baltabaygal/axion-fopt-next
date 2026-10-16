# Fig.2 physical-parameter rerun

H*/m=0.05,beta/H*=4,v=0.15,theta_i=2,L*m=8*pi,N=128,sharp
arrival mass (mf=0,mt=f=omega=1). Seed 0,100000 published-method tracers;
fixed outer nucleation dtau=0.01,indexed engulfment, births tau=0..55.
This matches stated physical/figure parameters, not the unpublished original
seed or exact adaptive integrator. Do not force full conversion from tracer F.

Start physical tau=10,t=2.5 before the first birth tau=15.58. The field is
exactly constant with zero energy before that birth, so a finite numerical
start reproduces its earlier massless history. To respect the native a_initial=1
convention, use a gauge rescaling c=a_physical(t_start)=0.5:

    x_program=c*x_physical, tau_program=c*tau_physical,
    a_program=a_physical/c, u=c*(tau_physical-10), H0=0.2.

Program box is 4*pi but physical comoving box in paper's a*=1 convention is
8*pi. Program first-arrival map and times are scaled by c. Physical densities
are invariant under this coordinate gauge change; plotted t=H*tau_physical²/2.
Evolve to physical tau=63.25,t=100.014, du=0.005 (physical dtau=0.01).
Also run N64 at du=0.005,0.0025 for spatial/timestep energy diagnostics.

Project the large catalogue using cell-bounded exact first arrivals. Preserve
the original brute-force oracle and compare N8 maps and random N128 cells
within 2e-12 absolute time before field evolution. Validate normalized scale
and snapshots/energies against independent Python formulas at endpoints and
checkpoint; max local mass/arrival mismatch 0, energy difference <1e-10.
N64 halved-timestep final total energy difference must be <1%; spatial N64/N128
difference is reported, not forced to pass. Sharp switches need independent
temporal checks; a single realization does not establish statistical convergence.
No-FOPT comparison uses homogeneous nonlinear cosine RD ODE initialized with
regular early theta_i=2. Zero FOPT energy before the first birth is exact for
this catalogue, not fabricated missing data. Save catalogue, rate, projection
checks/times, native energies, physical parameters, seed/source/build identity.
