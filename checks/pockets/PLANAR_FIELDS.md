# Native one-dimensional coupled fields

2026-10-18. This is an actual CosmoLattice NDim=1 model using its PV2
integrator and the same scalar potential/expression wrapper as PocketModel.
It is not a prescribed-wall evolution, particle simulation or separate
home-built field integrator. Both phi and a evolve and backreact.

The new `line_model.h` and `line_evolve.cpp` leave the existing 3D driver
unchanged. A CMake target `pocket_line_evolve` is available. The lightweight
validation executable is `build/pockets-line-light-01/pocket_line_evolve`;
that directory's build.json records its source hashes, compiler/link commands
and executable hash. A single unoptimized translation unit reused existing
CPU libraries and compiled in 1.76 seconds; no upstream files were edited.

## Domain and units

A periodic line, canonical kinetic terms, Minkowski cosmic time, no expansion,
PV2 only. Its energy is per transverse area. The evolution has one spatial
Laplacian; there is no spherical 2/r term. This cannot determine 3D pocket
stability, angular stress isotropy, spherical mass-radius scaling or bubble
collision efficiency. It is intended for local scattering and wall response,
to be compared against the independent continuum ODE reference in
FORMATION_SCATTERING.md.

Initial input is N rows of phi,a,phi_dot,a_dot at cell centres. Inputs must
be finite, without extra trailing values. N is limited to 4..16384; dt/dx<=.25,
positive integral timesteps, static background, and PS_type=2 are required.
The framework requires PS_type=2 for its native non-3D mode; this is a
configuration requirement, not a modified EoM. No spectrum estimator is
used here. Potential parameters are lambda, delta_v, mass_fv, mass_kappa,
axion_fa. The model derives mTV=mFV exp(kappa/2).

Outputs are initial.txt/final.txt, actual parameters and history.csv. The
history gives time, total/phi/axion energies, signed Q_a integral over the
line, and dt times the full mixed-Hessian lattice frequency bound. Each
potential-dependent diagnostic uses the formation kernel. Periodic global
currents cancel, so dE_a/dt=integral Q_a; dE_phi/dt=-integral Q_a. Future
spatial reflection/transmission measurements need explicit local windows or
matched control subtraction, not just these global energies.

## Validation evidence

`validate_line_backend.py` requires a matching frozen build manifest and
checks the executable/source hashes before running. Evidence:
`outputs/pockets/formation-line-validation-01`, including the declared gates,
raw initial states, histories, final fields, independent final reference,
parameter files and executable provenance. All runs use one thread.

N=32, L=8, dx=.25, t=0..0.2, dt=.01/.005/.0025. Smooth periodic
inhomogeneous fields and both momenta activate gradients, finite false mass,
nonlinear exponential forces and nonzero exchange. An independent DOP853
integration uses the same 1D nearest-neighbour stencil. Final full-field
absolute errors are 6.2793e-6,1.5697e-6,3.9243e-7, showing second-order time
convergence. The finest relative energy drift is 9.9449e-8; the integrated
axion exchange residual is 9.7324e-8. Snapshot energy agrees within 1.92e-16
of initial total energy. Initial transfer and an integer periodic translation
match exactly. These checks completed in under one second.

This verifies a low-cost route to coupled local scattering. It does not
establish spatial resolution or a scattering result yet. Next prepare a
well-separated periodic kink/antikink pair and a spectrally documented weak
axion packet. Check the empty-wall control, finite bandwidth, initial wall
overlap, causality/image times, time/space refinement and low-amplitude
limit against the continuum reference. Moving walls and finite-amplitude
backreaction follow as separately controlled cases. Keep all runs small
under the current laptop constraint.


Spectrally documented initial data and matched empty-wall inputs are now
prepared, with no solver launch. See
[PLANAR_SCATTERING_PREPARATION.md](PLANAR_SCATTERING_PREPARATION.md) for
localization, dispersion contamination, boost signs and image-time limits.


Short stationary/moving empty-wall controls now pass at two grids in 2.30
seconds total solver time; see PLANAR_EMPTY_CONTROL.md. An independent
moving-window energy diagnostic also passes instantaneous continuity tests.
Full scattering-time controls and compiled local diagnostic outputs remain
pending; no reflection/transmission fraction has yet been measured.


Optional passive local window outputs are now source-integrated and validated
in actual compiled evolution: see PLANAR_WINDOW.md and
`outputs/pockets/planar-window-validation-01`. Enabled diagnostics leave
fields and existing history values exactly unchanged in the small checks.
