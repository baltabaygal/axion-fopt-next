# Field mass-function comparison: phi^2 and phi^4

Declared before changing or evolving the model, 2026-10-17.

The phi^2 studies show early shrinkage and substantial axion-to-dilaton
energy transfer, surviving history refinement and mean-wall preparation.
Test the role of the nonzero false-vacuum mass curvature by adding the
explicit alternative

    V(phi,a) = U(phi) + 1/2 m_TV^2 phi^p a^2, p in {2,4}.

This is a toy coupling family, not the exponential instanton mass of the
paper. In p=4 the coupling's first and second phi derivatives vanish at
phi=0, unlike p=2 whose second derivative is m_TV^2 a^2. Both retain
m_a(0)=0 and m_a(1)=m_TV. Turning surfaces, wall forces, spectra and state
counts also change; a comparison cannot be attributed solely to trapping.
No lifetime fitting or unchanged-paper-physics claim is permitted.

## Implementation and checks first

Preserve the default p=2 potential and force and archived baseline build.
Validate p=4 independently by finite-differencing the Python potential,
compiled-versus-DOP853 field evolution, time refinement and energy-accounting
checks including coupling and far field. Use the full largest eigenvalue of
the local two-field Hessian plus the lattice Laplacian bound for the timestep
guard. The old maximum-diagonal estimate omits the mixed curvature.

At phi=0 the p=4 coupling curvature must be zero for every axion amplitude;
at phi=1 its axion mass squared remains m_TV^2. Reject exponents other than
2 or 4, including fractional input. All source/build identities are saved.

## Pilot comparison after validation

N=128, L=64, dx=.5, dt=.0125, t_end=80, history every .5,
slices every 4, full initial/final states. R0=20, w=2, lambda=2,
DeltaV=.01, m_TV=2. Same tanh initial wall, no wall velocity.
Each coupling uses its own finite-mass bound modes with .15<=omega<=.6,
radial h=.01, r_max=31, r_cut=30. Gaussian equipartition field and momentum
coefficients as in the original baseline, seed=0, with initial lattice
axion energy 2680.825731. Amplitudes and mode counts must be recorded.
This changes the gas shape and is an explicitly recorded confound.

Numerical gates: exact transfer, completion, energy drift<5e-4, record full
slice contour diagnostics without dropping missing brackets or thin patches.
An unresolved pilot is a screening result, not a trusted lifetime.
Compare radius loss, axion depletion, far axion/dilaton storage and spectra.
Require at least a factor-two reduction in axion depletion and radius loss
to identify a useful improvement. Report results either way. Subsequent
refinement and fixed-state box checks are required before long-time runs.

Default p=2 must reproduce the prior balanced trajectory within radius .001
and far stored fraction 1e-5 at all shared times; the stronger guard can
stop a run if necessary, which is reported rather than relaxed.
