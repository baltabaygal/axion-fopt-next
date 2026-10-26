# Fixed-sphere energy-current monitor

Declared before compiled implementation, 2026-10-17. A long pocket lifetime
needs energy-flow evidence, not just far-field storage in a periodic box.
Add an optional diagnostic without changing field evolution or boundaries.

For static canonical fields on the existing nearest-neighbour lattice,
split each bond's gradient energy equally between its endpoints. For either
field q with canonical momentum p, the outgoing face current is

    J_(i+1/2) = -(q_(i+1)-q_i) (p_(i+1)+p_i)/(2 dx).

Sum J*(inside_i-inside_(i+1))*dx^2 over all positive-axis faces to get
the net outward power from a fixed, box-centred sphere. This identity
follows by differentiating the symmetric bond energy and using the
semi-discrete field equation. It retains incoming power, including any
periodic return. Do not call its positive or net integral irreversible loss.

Allocate the coupling potential to the axion sector. Its signed exchange is

    Q = integral_inside (p_mass/2) m_TV^2 phi^(p_mass-1) a^2 pi_phi dV,

positive into axions. The exact semi-discrete identities are

    dE_phi/dt = -P_phi - Q,
    dE_a/dt   = -P_a + Q.

The old region energies use forward-bond assignment. They are retained
unchanged. The new monitor energy uses symmetric bonds and is explicitly
separate. All new columns appear only if 0<monitor_radius<L/2. The default
monitor_radius=0 leaves the existing output columns and dynamics unchanged.

## Numerical gates before physical use

An independent Python implementation must obey the identities by
directional finite differences for both phi^2 and phi^4 coupling.
Compiled monitor energies, powers and exchange at initial/final states
must match Python to relative 1e-12 (fixed total-energy scale).
Integrated continuity is tested at dt=.005 and .0025, t=1, N=6, L=6,
radius=2.2, with every-step histories; the total and sector residuals must
decrease at least by a factor 3.5 and be below 1e-4 at the coarse step.
Reject negative, nonfinite and >=L/2 radii before making native outputs.

Future physical runs must refine their output cadence separately for flux
integrals. A validated instantaneous monitor does not itself establish a
converged small loss rate, isolation, or long pocket survival.
