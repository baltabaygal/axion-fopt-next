# Fixed-action field preparation comparison

Declared before either compiled-field run, 2026-10-17.

The balanced Gaussian pilot passes its early history-refinement and box
comparisons, but transfers substantial axion energy to the dilaton and has
unresolved patches. Test whether a self-consistent mean wall helps before
changing the field potential. This is a new classical initial-state family.

## Independent preparation

`relax_reference.py` minimises the spherical time-averaged energy

    E[phi] = 4 pi integral r^2 [phi_r^2/2 + U(phi)] dr
             + sum_(l,n) (2l+1) I_ln omega_ln[phi],

with the same phi^2 axion mass and scalar potential as the compiled solver.
Eigenfunctions obey -u''+[l(l+1)/r^2+m_TV^2 phi^2]u=omega^2 u.
The Hellmann-Feynman derivative supplies the mean force. Directions are
isotropic only at the ensemble level. The thermal/equipartition choice
I_ln=(2680.825731/93)/omega_ln(initial tanh) defines conserved actions for
this preparation, not the actions of the earlier Gaussian realisation.

The h=0.1 and h=0.05 references have R_half=19.65857 and 19.65880,
energies 3860.664 and 3860.821, and force-density residuals <8e-14.
These are reference equilibria, not field evolution or stability results.

## Compiled field pair

Use h=0.05, r_max=40 reference, and the same actions in both cases:

- `tanh`: phi=[1+tanh((r-20)/2)]/2, with its own bound eigenfunctions.
- `relaxed`: minimised phi with its own bound eigenfunctions.

For each real spherical harmonic, choose a phase uniformly on [0,2pi)
using seed=0 in (l,n,m) order, and set

    q_lnm=sqrt(2 I_ln/omega_ln) cos(phase),
    p_lnm=sqrt(2 I_ln omega_ln) sin(phase).

Each harmonic then has the specified action and continuum energy I*omega.
Fix eigenfunction sign by its maximum absolute sample being positive.
This is a coherent, high-occupation classical field construction, not a
Bose thermal gas. Degenerate harmonic interference remains in a single
realisation. Do not claim exact local isotropy or thermalisation.

No lattice energy renormalisation. Sample each continuum state at N=128,
L=64, dx=0.5, dt=0.0125, to t=80; history every 0.5, slices every 4.
Axion samples at r>=30 are zero (exponentially small tail); phi=1 beyond
r_max. lambda=2, DeltaV=.01, m_TV=2; original mode labels selected from
0.15<=omega<=0.6, recomputed in each profile. Full initial/final states.

Numerical gates: exact transfer, completion, energy drift<5e-4. Diagnose
full slice contours, reporting all missing brackets and widths<4 cells;
these cannot be excused by an energy gate. No continuum lifetime claim.

Report R(t), axion-to-dilaton energy transfer, far stored axion energy,
and radius/energy loss relative to each initial state. A useful improvement
requires at least halving both radius loss and axion energy depletion at
t=80, with no failed numerical gate. Otherwise the mean-profile adjustment
alone is insufficient. Even improvement requires later grid, box, phase
and perturbation checks before the ten-crossing-time goal is established.
