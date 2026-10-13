# Controlled field-evolution checks

This workspace investigates simple limits before more complicated FOPT runs:
homogeneous fields, then a single bubble nucleating at the center
(L/2,L/2,L/2) of a periodic cubic volume L^3. Non-expanding comparisons are a
possible next case. Implemented cases are the homogeneous radiation-era checks
below and the [non-expanding centered bubble](single_bubble/README.md).
The [stochastic nucleation workspace](stochastic/README.md) defines the next
geometry/statistics investigations and their validation order.

## Homogeneous cosine potential and the Bessel reference

We evolve V/f^2=m^2(1-cos(theta)) with constant mass everywhere and no bubbles.
The radiation background is a(tau)=H_* tau, t=H_* tau^2/2. In cosmic time:

    theta_tt + 3/(2t) theta_t + m^2 sin(theta) = 0.

Appendix B, Eqs. (B1)-(B2), of
[Baltabay, D'Eramo and Vaskonen (2026)](https://arxiv.org/pdf/2607.01333)
gives the regular small-angle solution after replacing sin(theta) by theta:

    theta_B(t) = theta_i 2^(1/4) Gamma(5/4) (mt)^(-1/4) J_(1/4)(mt).
    dtheta_B/dt = -theta_i 2^(1/4) Gamma(5/4) m (mt)^(-1/4) J_(5/4)(mt).

This is a linear reference, not an exact solution of the cosine potential.
An independent high-accuracy DOP853 integration of the full cosine equation
separates lattice discretization error from the small nonlinear correction.
The linear ODE is also integrated independently to check the Bessel formula
and its time/derivative mapping.

The adapter requires a_initial=1 and a finite conformal start. We choose
H_*=m=omega_star=f_axion=1, tau_origin=1, so t_initial=0.5,
u=tau-1, a=1+u and t=(1+u)^2/2. Both theta_initial and
theta_prime_initial=a_initial*dtheta_B/dt are evaluated at t_initial.
theta_i denotes the regular solution's extrapolated t->0 amplitude, not the
field value at this finite start. We do not impose a zero derivative at t=0.5.

The fixed physical box is L=2*pi/m. The grid sweep Ngrid=4,8,16,32 uses
theta_i=0.01, dt=0.002 in conformal/program time and u_final=8
(mt_final=40.5, about six oscillations). Temporal refinements at Ngrid=8
use dt=0.001 and 0.0005. A theta_i=0.001 run at the finest timestep checks
the small-angle limit. No stochastic sampling is involved: every arrival-map
cell is already converted at the start.

## Run from the repository root

The cached adapter is used by default; build instructions are in
[docs/cosmolattice.md](../docs/cosmolattice.md). Requires the project plus the
validation and plots extras (NumPy, SciPy and Matplotlib).

```sh
PYTHONPATH=src MPLCONFIGDIR=/tmp/afopt-mpl python3 checks/homogeneous.py \
  --output outputs/homogeneous-checks
```

The output directory must be fresh. Optional `--executable` selects another
compiled adapter. Full run inputs, logs, initial/final cell snapshots,
histories, reference curves, figures, and summary.json are saved there.
The script writes all measured gates to summary.json and exits nonzero on a
failed gate. It never uses Python assertions for acceptance.

Errors are max absolute differences divided by the fixed theta_i, avoiding
singular relative errors at field zero crossings. Predeclared gates:

| Check | Gate |
|---|---|
| Independent linear ODE versus Bessel | normalized field and momentum error < 1e-8 |
| Cosine lattice versus cosine ODE, dt=0.002 | normalized field error < 3e-4 |
| Finest timestep versus cosine ODE | normalized field error < 2e-5 |
| Small-angle baseline versus Bessel | normalized field error < 5e-4 |
| Grid-to-grid field/momentum histories | normalized difference < 1e-10 |
| Initial/final local field/momentum peak-to-peak spread | normalized spread < 1e-10 |
| Temporal refinement | coarse/fine field-error ratio > 10 for fourfold dt reduction |
| Nonlinear correction | decreases by > 80 for tenfold amplitude reduction |

Momentum means pi=a^2*dtheta/du. The reference momentum uses the same
canonical definition. Scale-factor, time, mass and initial-state transfers are
also checked. Mean-square histories provide a roundoff-limited diagnostic of
spatial spread during the run; exact local spread is measured at the saved
initial and final states.

Plots use time-series overlays and normalized residuals for analytic agreement,
and log-log error versus timestep for refinement. Raw reference CSVs retain
the entire sampled time series. These are finite-duration numerical checks;
grid independence here does not establish spatial convergence for a bubble.

See [results/README.md](results/README.md) for the recorded run once generated.

## Runtime optimization

[benchmark_runtime.py](benchmark_runtime.py) compares saved baseline and
optimized executables on the same N=128, v=0.5 bubble case: 1280 steps,
unchanged double precision, wall profile, diagnostics every two steps and two
binary snapshots. It alternates run order and checks every final cell and
history column against the original bubble output. The executable wall time
includes initialization and output; Python comparison and FFT analysis are
outside the timed region. macOS peak RSS is captured separately.

```sh
PYTHONPATH=src python3 checks/benchmark_runtime.py \
  --baseline outputs/runtime-128-baseline/profiled_evolve \
  --optimized build/cosmolattice/axion_fopt_evolve \
  --output outputs/runtime-128-cached
```

The default repeats are three alternating pairs. Use a fresh output directory.
The historical baseline binary is intentionally kept in ignored outputs/;
it must be preserved or rebuilt from the archived pre-optimization source
before repeating this comparison. See [runtime_results/README.md](runtime_results/README.md)
for recorded timings, provenance and numerical differences.
