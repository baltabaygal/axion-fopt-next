# First F(t) comparison: frozen before executing ensembles

Primary observable: geometric false-vacuum fraction versus analytic survival,
with independent probe points and independent catalogue realizations.
Also retain the paper method's internal active/total tracer estimator.

Two geometry-only cases, c=1, periodic cube, common constant wall speed:

| Case | q(tau) | L | v | tau_start | tau_end |
|---|---|---|---|---|---|
| Static control | 4 | 3 | 0.5 | 0 | 1 |
| Radiation era | Gamma(t) a^4; H*=beta=1, a=tau, t=tau^2/2, t*=0.5 | 3 | 0.5 | 0 | 2.5 |

The second case is a controlled slow-rate diagnostic, not a full reproduction
of the paper's Fig. 1 parameter scan. All radii remain strictly below L/2.
The physical radiation background starts at cosmic time zero; the formal
negative-time lower limit in the paper is not extended through the Big Bang.
No nucleation cutoff at t* is inserted.

Reference: F=exp[-(4*pi/3) v^3 integral_0^tau q(s)(tau-s)^3 ds].
Independently integrate the same survival in cosmic time using Gamma*a_birth^3
and the comoving radius. Compare the physical q reference with the exact
piecewise-linear interpolant reference used by the new sampler.

Initial ensemble: 24 seeds (0..23), 2048 independent uniform probes per
realization split into two halves, independently seeded from nucleation RNG.
Poisson rate grids: 65 and 257 nodes. Tracer population: 512 and 2048,
with a finer-time 2048-tracer comparison. Static tracer outer dt=0.02,0.01;
radiation outer dt=0.05,0.025. Output times coincide with both timestep grids.
These modest tracer populations deliberately expose finite-population effects;
neither substitutes for the paper's 100,000-tracer runs.

Primary acceptance gates, declared before results:

- Conformal and cosmic analytic F curves differ by <1e-9.
- Fine Poisson interpolant F differs from physical-rate F by <0.002.
- Poisson ensemble mean F agrees with its interpolant reference at every
  sampled time within 4 ensemble standard errors +0.005 absolute F.
- Two probe-half ensemble means differ within 4 paired standard errors +0.005.
- All geometric and internal survival curves are bounded and nonincreasing.
- For identical candidate events, rejected and all-candidate catalogues have
  matching arrival maps (absolute tolerance 1e-12); accepted births are causal.

Tracer/reference differences, tracer-population dependence, timestep dependence
and number counts are diagnostics, not forced-to-pass gates. Do not alter
parameters or relax thresholds after seeing an unfavorable outcome. The
ensemble uncertainty is uncertainty in the mean, not the realization spread.

The source is pinned to GitHub commit 1b67fec and blob ace7473. Load only the
paper's bubble/tracer classes and geometry kernels, with JIT disk caching
disabled. The geometry-only driver reproduces its nucleation ordering while
replacing the field-dependent outer timestep with the specified fixed dt.
Save the original source under ignored outputs/ alongside its hash. Preserve
the published source and do not import or execute its field engine.
