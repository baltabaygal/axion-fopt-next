# Instructions for the complex-PQ string and wall workspace

Read README.md here and ../AGENTS.md before changing this workspace. This
directory owns the modest-grid post-inflationary PQ experiment requested by
the user. It is separate from prescribed FOPT bubble simulations.

- Implement the complex field as two canonically normalized real components.
  Preserve the radial mode and derive both forces from the documented potential.
  Do not replace dynamical strings with arrival maps or an angle-only field.
- Follow arXiv:2609.20801v1, Eqs. (1)-(4), (45), and Sec. IV. State every
  deliberate simplification; reduced volume alone does not preserve the paper's
  hierarchy of horizon, correlation, core and wall scales.
- Use the actual compiled CosmoLattice/TempLat backend for reported lattice
  evolution. Keep the existing one-scalar adapter and its tests working. Put
  experiment-specific model/build/driver code here where practical; generated
  builds belong under the ignored root build/ directory. Never edit upstream
  vendor source. PV2 is the initial backend choice and differs from the paper's
  leapfrog implementation; validate it independently.
- Use conformal time and comoving coordinates. Document the scale-factor
  normalization, canonical momenta, lambda, f_a, susceptibility law and every
  time/stage conversion. Evaluate prescribed coefficients at actual force stages.
- Treat r=|Phi|=0 explicitly. Derive a finite Cartesian force for the paper's
  r^2*(1-cos(theta)) tilt, including its limiting value at the origin. Do not
  silently add an epsilon, substitute a linear tilt, or differentiate atan2.
  Current/spectrum observables need their own documented treatment of cores.
- Start with analytic/ODE and deterministic periodic defect checks before a
  random network. A periodic string pair or loop is appropriate; a lone straight
  net-winding string is not a neutral periodic initial condition.
- For random initial conditions record the covariance/filter convention,
  variance definition, correlation length, seed and RNG. Cooling and core-growth
  equations/schedules must be explicit. Verify referenced preparation details
  before claiming the paper's initialization has been reproduced.
- Check dt/dx and the largest a*m for both radial and angular modes. Report
  core and wall widths in cells throughout the run. Stop quantitative claims
  when defects are unresolved; do not infer stability from CFL alone.
- Predeclare numerical gates before executing each comparison. Use an independent
  small-grid semi-discrete ODE reference, temporal refinement, and separate spatial
  refinement with fixed physical parameters. Freeze initial data appropriately
  when comparing discretizations; reseeding alone does not match different grids.
- Measure topology with periodic wrapped phase increments and verify known
  winding configurations. Distinguish raw plaquette/link counts from calibrated
  string-length/wall-area estimators and from continuum convergence.
- Record energy components and the explicit source/expansion terms. Ordinary
  energy conservation applies only to static-background, time-independent tests.
  Cooling and core growth are preparation stages, not physical production data.
- Save parameters, exact initial states, seeds, source/build identity, histories,
  diagnostic definitions and machine-readable checks. Raw outputs go into fresh
  directories under outputs/pq_strings/; compact evidence may live in results/
  here. Never overwrite earlier runs.
- Keep initial runs modest: N=32 pilot, then N=64 only after correctness checks.
  Larger runs, distributed execution and paper mass/abundance predictions are
  separate work. Do not claim scaling from a small box or one realization.

The user authorizes this experiment's routine implementation and reversible
checks. No extra approval step or agent delegation is required.

## Authorized capped-mass toy phase

The user explicitly selected a modified small-grid toy after the fixed-grid
feasibility audit. Keep the full PQ fields and potential, but allow the documented
smooth capped susceptibility. Label this phase as finite-box collapse/yield work;
never claim paper reproduction or infer a physical axion mass from its yield.
Keep the default paper-inspired power law available and tested. A late homogeneous
mode completion is separate from the current-only estimator and needs explicit
small-angle and radial-vacuum gates. A 128³ matched spatial refinement is authorized
for this toy after the validated 32³/64³ stages; no larger grid is authorized here.
