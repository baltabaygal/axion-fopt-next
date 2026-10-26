# Fourth-order spatial reference: lower-cost accuracy path

2026-10-19. No solver was launched and no compiled evolution source or
vendor file changed. The current actual line backend still uses its original
second-order nearest-neighbour stencil. The moving fixed-speed 1% waveform
gate remains failed at 5.568% on the finest recorded grid.

This checkpoint develops and checks a consistent alternative before
implementing it. Independent formulas live in `planar_stencils.py`; raw
audit evidence and hashes are in
`outputs/pockets/planar-fourth-order-reference-01/report.json`.

## Spatial operator and matching energy

For hop coefficients c1=4/3, c2=-1/12,

L4 f_i = sum_(l=1,2) c_l (f_(i+l)+f_(i-l)-2 f_i)/dx^2.

The endpoint-shared gradient energy density is

rho_grad,i = sum_l c_l [(f_(i+l)-f_i)^2+(f_i-f_(i-l))^2]/(4 dx^2).

The summed gradient energy is nonnegative. Its Fourier symbol is
4 sin^2(k dx/2) [1+sin^2(k dx/2)/3]/dx^2, ranging from zero to
16/(3 dx^2). Individual long-bond contributions and local split densities
need not be positive. They must not be interpreted as separate physical
negative-energy modes or used uncritically for energy-weighted geometry.

The Hamiltonian derivative gives exactly -L4. Replacing only the force
operator while keeping nearest-neighbour energies or currents would violate
accounting. The maximum kinetic frequency contribution also changes from
4/dx^2 to 16/(3 dx^2), so the timestep guard must change consistently.

## Long-bond currents

For each hop l, the edge current is

J_(i,i+l) = -c_l (f_(i+l)-f_i)(pi_i+pi_(i+l))/(2 dx).

For an arbitrary periodic diagnostic window W, outward lab flux is
sum_(i,l) J_(i,i+l)(W_i-W_(i+l)). Boundary advection is the same weighted
rho*dW/dt integral as before; corrected flux is lab flux minus advection.
The axion source remains kappa V_ax pi_phi. The signed second-neighbour
edges are included explicitly rather than folded into the old face current.

Three targeted tests pass: Fourier symbol/positivity and continuum order;
independent finite-difference Hamiltonian derivative versus force; and
instantaneous total/phi/axion moving-window continuity with nonzero exchange.
The order-2 branch also matches the established independent diagnostic.
These are mathematical/reference checks, not compiled adapter tests.

## Measured reference benefit

Sample the same continuum degenerate kink pairs as the existing controls,
with fixed lambda=8, width=1, L=512 and speeds 0/.3. Check
(1-u^2)*L(phi)-U_phi; a continuum boosted kink obeys zero residual.
No field evolution or speed/phase fitting is performed.

| Moving u=.3 grid | dx | Order-2 residual | Order-4 residual |
|---:|---:|---:|---:|
| 2048 | .25 | .011326 | .0011515 |
| 4096 | .125 | .0028837 | .000082696 |
| 8192 | .0625 | .00072954 | .0000052860 |

At the current middle grid the residual improves by ~34.9 times at fixed
spacing. The stationary residual improves by ~38.4 times. Separate declared
ratio gates confirm second-order and fourth-order behaviour for both speeds.
Linear incoming and reflected-mode dispersion estimates are also saved in
the audit. This motivates the approach, but a residual improvement does not
prove that a moving packet's accumulated waveform error will pass 1%.

## Required compiled integration before a physical claim

Add an explicit optional planar Laplacian hook in a user-owned adapter
extension, preserving the default order-2 behaviour for other models and
archived runs. Do not change vendor source or disguise gradient corrections
as a potential derivative. The hop-2 operator needs two ghost cells rather
than the existing one. Its energy, passive-window currents and frequency
bound must be wired together.

Then check compiled forces against the independent semi-discrete reference,
energy/current identities, time and spatial order, and default-order numerical
regression. Only after those checks run the existing moving benchmark with
its same physical data, fixed-speed reference and unchanged 1% gate. No
initial velocity calibration or fitted phase may be used to conceal the
original error. Large grids and 3D production runs remain deferred under
the laptop constraint.
