# Field-frequency hierarchy pilot

Contract: [FREQUENCY_EXPERIMENT.md](../FREQUENCY_EXPERIMENT.md), declared
before the runs. Raw fields are in `outputs/pockets/frequency-pilot-01`;
the authoritative reconstructed analysis is in
`outputs/pockets/frequency-analysis-01/report.json`.

All three compiled field evolutions completed t=80. The first runner's
JSON-writing stage failed on a NumPy boolean; its tensor reshape for the
asphericity was also corrected before results were used. Reports were
reconstructed from saved initial/final states and histories into a fresh
directory, without repeating or changing any evolution. Original source
snapshots remain alongside the raw suite. The corrected covariance assembly
is checked against a rotated ellipsoid and an isotropic tensor.

## Parameters and result

Static a=1, N=128, L=64, dx=.5, dt=.0125, lambda=2, DeltaV=.01,
m_TV=2, phi^4 coupling, R0=20, w=2, seed=0. Each case has initial lattice
axion energy 2680.825731, its own finite-mass bound modes, radial h=.01,
r_max=31, r_cut=30. The Gaussian classical high-occupation field has not
been thermalised. All include the fundamental mode via omega_min=.1.

| Case | omega_max | Modes | Final R | Axion-sector depletion | Far axion storage / E_a(0) | Far dilaton storage |
|---|---:|---:|---:|---:|---:|---:|
| soft | .35 | 17 | 18.517 | 13.159% | 2.497% | 198.02 |
| intermediate | .45 | 46 | 17.828 | 9.366% | 2.069% | 181.29 |
| broad | .60 | 99 | 18.123 | 11.506% | 1.840% | 206.39 |

All exact-transfer, completion, initial energy-accounting and global
energy gates pass. Maximum relative energy drift is <7.3e-6.
No case passes the predeclared <=5% depletion candidate threshold.
Every case also fails the wall diagnostic due to unresolved or ambiguous
points. Final full-3D minimum measured thicknesses are 3.75, 4.01 and
4.29 cells respectively, with 118, 316 and 423 ambiguous brackets.
No missing bracket is dropped to obtain a pass.

The softest case reaches asphericity 1.029, compared with .813 and .527
for intermediate and broad. Changing the band also changes mode count,
amplitude and pressure fluctuations. There is no monotonic decrease in
depletion as the initial pair-sum ceiling falls below the exterior dilaton
mass 1.0583. The data therefore do not isolate a quadratic-source threshold
as the sole cause of radiation, nor establish an improved long-lived state.

The characteristic gas frequency and complete diagnostics are in the
machine-readable report. Far energies here are storage, not cumulative
escape or radiation flux. Periodic return and incomplete wall resolution
still limit interpretation.

## Next experiment

Increase the dilaton mass itself with a narrower, resolved wall while
retaining a populated broad axion band. The
[stiff-wall contract](../STIFF_WALL_EXPERIMENT.md) records the changed
tension and energy normalisation, paired grids, and the fixed-sphere
energy-current diagnostic. A successful screening result still needs the
unchanged t>=400 survival, perturbation, convergence and isolation gates.
