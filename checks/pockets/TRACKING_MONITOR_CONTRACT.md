# Energy current around a moving field pocket

Declared before implementation, 2026-10-17. The stiff-wall pocket moves
out of a fixed sphere, so fixed-sphere export is not whole-pocket loss.
Develop a moving control-volume reference before extending lifetime claims.

Use a smooth phase weight eta(phi)=1-S((phi-.1)/.8), where the argument
is clamped to [0,1] and S(s)=6s^5-15s^4+10s^3. Compute a circular
centroid to choose an unwrap cut, then the arithmetic centroid of this
localised weighted pocket in unwrapped coordinates. Reject ambiguous
centroids, rather than silently picking a position for a uniform/multiple
phase distribution. In the static model:

    eta_dot = eta'(phi) pi_phi,
    c_dot_i = sum eta_dot (x_unwrapped_i-c_i) / sum eta.

For a compact radial window W about c, W=1 inside R_in and W=0 outside
R_out<L/2, use the same quintic ramp. W_dot=-W_r (x-c).c_dot/r.
Symmetric bond energies rho and the existing exact lattice face current J
give the exact semi-discrete identity

    d sum(W rho) dV / dt = -sum_faces J (W_i-W_j) dx^2
                           +sum_cells W_dot rho dV.

Thus the net current out of the moving window is lab-frame current minus
the advective term. In each field sector add the existing signed coupling
exchange Q; it cancels in the total. Do not call fixed or moving net
current irreversible escape, since incoming waves remain possible.

The reference uses R_in=2.5,R_out=3.5 on a small L=8 test box; later
pocket runs would use R_in=26,R_out=30 on L=64 only if the bound material
lies within the inner window. Output both lab and advective terms so that
translation corrections remain inspectable.

## Gates before compiled physical use

- Independent directional finite differences of the weighted energies
  must match total and sector currents plus coupling exchange to 1e-7
  relative (fixed initial weighted-energy normalisation), including actual
  centroid motion, for p=2/4 and vacuum_stiffening=0/12.
- Verify invariance under integer-cell periodic translations, including a
  pocket straddling the box edge, to 1e-12 relative.
- A compiled implementation must independently match energies, currents,
  centroid and velocity against this reference, and exhibit second-order
  integrated-continuity/time refinement before physical use.
- The existing fixed monitor and field dynamics must remain unchanged by
  default. A finite moving-window storage estimate alone does not prove
  isolation, rest mass, stability or a long lifetime.
