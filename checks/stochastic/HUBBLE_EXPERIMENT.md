# Hubble-rate validation (declared before execution)

Geometry only, radiation domination: a=H*tau, t=H*tau^2/2,
Gamma=H^4 exp[beta(t-1/(2H))], q=Gamma*a^4. Start at tau=t=0.

Two comparisons: H=0.5,1,2 with beta/H=1 (scale covariance), and
H=0.75,1,1.5 with beta=1 (a physical change in beta/H).
The H=beta=1 ensemble is shared between comparisons.
Keep H*L=3, v=0.5, H*tau_end=2.5, H*dt=0.05; thus
H*t_end=3.125 and maximum radius is below L/2 for every case.
This is a Hubble-sized box comparison, not fixed physical box length.

Each case: seeds 0..23, 2048 independent probes, Poisson rate table with
257 nodes, and the pinned paper method with 2048 tracers. Probe RNG uses
SeedSequence([92842, seed]), rescaled by L, so normalized probe positions
match across H. Save catalogues, counts, rate tables, internal and geometric
F, analytic curves, ensemble standard errors and workflow wall times.

Acceptance: conformal/cosmic analytic agreement <1e-9; 257-node linear-rate
reference error <0.002 and its difference from 513 nodes <0.002;
Poisson mean within 4 ensemble standard errors +0.005 at each observation;
probe halves within 4 paired standard errors +0.005; bounded monotonic F;
accepted Poisson births causal; phantom arrival maps agree within 1e-10 in
H*tau. Fixed-beta/H analytic curves agree within 1e-10, and identical-seed
Poisson normalized arrival maps agree within 1e-10 across H.
Tracer differences remain diagnostics, not forced acceptance gates.

Timings include catalogue generation, independent geometry, validation and
compressed output, and the first tracer batch includes JIT compilation.
These are small ensemble checks, not isolated kernel benchmarks or coupled
field runs. Preserve previous experiment outputs.
