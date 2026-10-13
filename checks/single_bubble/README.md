# Centered periodic bubble — Lee et al. (2024)

Source: [arXiv:2402.09501v2](https://arxiv.org/pdf/2402.09501), Sec. IV F,
Eqs. (53)-(55), Table I. We implement its non-expanding harmonic equation,
not the cosine potential used in the homogeneous radiation-era check:

    phi_tt - laplacian(phi) + m^2(t,x) phi = 0,
    m^2 = m0^2 [1 - tanh((r-vt)/delta_w)] / 2.

One bubble starts at the center of a periodic cube. The false-vacuum mass is
zero. Set m0=phi0=1, phi(t=0)=1, phi_dot(t=0)=0, delta_w=0.15/m0,
L=2*pi/m0, and end at t=L/v. We run v=0.1,0.5,0.9 at Ngrid=128,
matching the smallest Table-I box. Larger published boxes are not run here.
The width is measured in the box frame. The periodic copies interact through
the scalar field; no separate wall-collision dynamics is introduced.

## Geometry and solver mapping

The center is (L/2,L/2,L/2), with cell-centered coordinates and periodic
nearest-image distance r. The arrival-map value is r/v. Setting
`wall_width=0.15` and `profile_wall_speed=v` makes the adapter evaluate
0.5*(1+tanh((t-arrival)*v/delta_w)), exactly the intended mass profile.
Its small nonzero initial core/tail is retained rather than replacing it with
an initially empty sharp bubble. The same nearest-center prescription handles
overlap with the periodic images, approaching a uniform true-vacuum mass.

Adapter settings: potential=quadratic, expansion=false, tau_origin=0,
omega_star=f_axion=1, so t=u=tau, a=1 and pi=phi_dot. PV2 evolves the
field with the existing periodic nearest-neighbor Laplacian. This matches the
equation and physical parameters; the paper does not specify its complete
time-integrator/stencil/start-grid implementation or supply raw outputs here.

The base step is L/v divided by ceil((L/v)/(0.2*L/128)). This reaches the
exact final time and satisfies dt/dx<=0.2. Each velocity's grid comparison
holds dt fixed. At v=0.5 we compare Ngrid=32,64,128 and halve dt at Ngrid=64
and 128. The finest grid has about 3.06 cells per wall-width parameter;
coarser runs are resolution diagnostics, not substitutes for the published grid.

## Run

From the repository root, after building the adapter:

```sh
cmake --build build/cosmolattice --target axion_fopt_evolve -j 4
PYTHONPATH=src MPLCONFIGDIR=/tmp/afopt-mpl python3 checks/single_bubble/run.py \
  --output outputs/lee-single-bubble
```

Use a fresh output directory. `--suite pilot` runs only the three N=32/64
v=0.5 cases; default `--suite paper` adds all N=128 cases. `--executable`
selects the compiled adapter. Individual input JSONs retain exact parameters;
logs, initial/final binary snapshots, histories, spectra and per-run observables
are kept in the output directory. Summary JSON is updated after each finished
run. A v=0.5, N=128 checkpoint is saved at the first step at or after mt=2,
before periodic wave overlap, to inspect initial depletion. Its actual time
is stored in the snapshot header.

## Observables and acceptance checks

Compute mean kinetic, forward-difference gradient and quadratic-potential
energies. The inside energy is the sum over r<vt divided by the entire box
volume. The mass profile injects energy: the correct balance is
E(t)-E(0)=integral(dt <0.5 phi^2 partial_t m^2>), rather than constant energy.
The source integral uses the saved time samples and trapezoidal quadrature.

The final number is computed with a full FFT normalized by Ngrid^3:

    n = sum_k (|pi_k|^2 + omega_k^2 |phi_k|^2) / (2 omega_k),
    omega_k = sqrt(k^2 + m0^2).

Include all modes and the zero mode. Do not double-count conjugate modes or
replace complex modulus squares with real-part squares. A second estimator uses
the stencil dispersion k_latt^2=sum_i [2 sin(k_i dx/2)/dx]^2. Its spectral
energy must match the real-space gradient energy by Parseval's theorem. The
continuum estimator is the one compared with the paper. Shell spectra report
dn/dln(k) and exclude k=0; the zero-mode number is reported separately.

Predeclared numerical gates:

- Relative Parseval energy residual < 1e-12.
- Relative energy/work balance residual < 0.005.
- Maximum final mass deficit < 1e-4.
- v=0.5 number change under timestep halving < 0.002.
- v=0.5 number change from N=64 to 128 < 0.02.

The independent DOP853 PDE test, Fourier zero-mode/plane-wave normalization
test, transfer checks and input validation are in tests/test_bubble_adapter.py.
Enable integration tests with AFOPT_COSMOLATTICE_EXECUTABLE:

```sh
PYTHONPATH=src AFOPT_COSMOLATTICE_EXECUTABLE="$PWD/build/cosmolattice/axion_fopt_evolve" \
  python3 -m unittest discover -s tests -v
```

The approximate piecewise velocity fit of Eqs. (56)-(57), originally measured
at L=16*pi, is shown for context. It is not an exact acceptance target for
L=2*pi. No raw-author-data or digitized-figure precision is claimed. A matched
smallest-box result does not establish the box-size scaling, the large-box
energy curve in Fig. 7, or the large-box spectrum in Fig. 10.

Recorded findings and figures are saved under results/ after the run.
