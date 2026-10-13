# Physics contract

Lengths are comoving and times conformal, with c=1. Gamma is a rate per
physical volume per cosmic time:

    dN = Gamma(t) a(t)^3 d^3x dt = Gamma(t(tau)) a(tau)^4 d^3x dtau.

The sampler accepts q(tau)=Gamma*a^4, interpolated linearly between supplied
nodes. It never multiplies q by P_FV before sampling: explicit rejection already
enforces exclusion. Resolve rapidly varying q by refining the input time grid.
There is no silent exponent clipping. The caller must supply finite intensities
and choose a time window whose candidate budget is tractable.

For radiation domination in the uploaded normalization:

    a(tau)=H_* tau; t(tau)=H_* tau^2/2; t_*=1/(2H_*).
    Gamma(t)=H_*^4 exp[beta(t-t_*)].

A constant physical speed obeys dx/dtau=v_w, hence
R_comoving=v_w*(tau-tau_birth). Arrival time is tau_birth+d_periodic/v_w.
Periodic nearest-image distances give the exact converted union on a cubic torus.
First-arrival maps can therefore be precomputed once for a fixed catalogue.

For a common constant speed, a rejected bubble born within an earlier bubble
cannot outrun its parent; the triangle inequality ensures phantom equivalence.
The unbounded-space expression is

    P_FV(tau)=exp[-I(tau)]
    I=(4*pi/3)*integral(q(s)*[v_w*(tau-s)]^3 ds).

The sphere-volume expression is valid on the periodic box while all radii are
below L/2. Beyond that use torus-ball volume, not 4*pi*R^3/3. Finite-box
correlations require ensemble errors, not independent-cell binomial errors.
Allowing bubble-specific speeds invalidates the phantom argument in general;
that extension requires a separate statistical derivation.

For dimensionless axion theta=phi/f:

    theta'' + 2*(a'/a)*theta' - laplacian(theta)
        + a^2*m^2(tau,x)*sin(theta) = 0.

The Python force API supplies m^2*sin(theta), without a^2. V/f^2 is
m^2*(1-cos(theta)). The solver supplies its own scale factor and unit conversion.
A changing prescribed mass injects/removes field energy, so ordinary closed
system energy conservation is not an appropriate acceptance criterion.
