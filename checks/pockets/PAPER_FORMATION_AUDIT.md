# Model and hierarchy audit for the accepted revised goal

Reference: [arXiv:2407.01676v2](https://arxiv.org/html/2407.01676v2).
Equation numbers below follow its PDF; the HTML numbers equations globally.
This is preparation, not a new coupled-field simulation.

## Model mapping

Let canonical displacement phi increase from FV=0 to TV=1, choosing the
dilaton excursion as the field unit. PDF Eq. (2.2) then motivates

    V(phi,a) = U(phi) + f_a^2 m_FV^2 exp(kappa phi)[1-cos(a/f_a)],
    m_TV/m_FV = exp(kappa/2).

The vacuum potential U is not fixed by the paper. An explicit benchmark
family can retain our tilted quartic, with its approximations reported.
This is not a claim to reconstruct a unique ultraviolet theory. Canonical
kinetic terms and the field excursion must be specified with any scale map.
The new independent `paper_reference.py` differentiates this potential and
its Hessian, including mixed terms. Tests check finite differences,
periodicity, vacuum masses, the quadratic limit and overflow rejection.
It is not yet implemented or validated in the compiled field backend.

Unlike the old phi^4 portal, exponential mass dependence produces dilaton
backreaction even at phi=0 for nonzero axion amplitude. A massless FV and
zero FV portal derivative are therefore not innocuous approximations to
this new benchmark. A finite-f_a cosine also provides intrinsic axion
nonlinearity absent from the old quadratic portal. Thermalisation must be
measured, not inferred from including a cosine.

## Moving-wall reflection and hierarchy

For an incoming mode with frequency omega, FV mass m, and normal momentum
k=sqrt(omega^2-m^2), a head-on wall approaching with speed u gives

    omega_wall = gamma(omega+u k),
    omega_reflected = gamma^2[omega(1+u^2)+2u k].

These follow Appendix A's Lorentz kinematics, without evolving particles.
Require the appropriate wall-frame frequency to lie below the exterior
threshold before treating reflection as efficient. Later encounters and
compression need their own spectral checks. With omega=.6,m=0,u=.999,
the wall-frame frequency is about 26.8: the old exterior mass 2 cannot
support the assumed reflective regime at that velocity. Lab-frame
omega<m_TV alone is insufficient for a contracting-pocket experiment.

For our quartic family's degenerate reference sigma=1/(3w), Rc=3sigma/dV.
The completed lambda=32,w=.5,dV=.01,R0=20 toy has Rc/R0=10, versus .01
in Fig. 4. It started at rest and predominantly tested a different regime.
Rc here is a thin-wall reference, not a measured bounce radius.

For a boosted tanh initial wall its 10--90 thickness estimate is
2 atanh(.8) w/gamma. Thus a uniform cube needs at least

    N >= (L/R0)(R0/w) 4 gamma/[2 atanh(.8)]

to put four cells across that initial estimate. Combined thin-wall,
large-radius and highly relativistic limits are beyond this 16-GB host's
uniform 3D capacity. `audit_paper_hierarchy.py` makes those assumptions and
costs explicit; it does not turn an inaccessible hierarchy into a passed
experiment. Next select an accessible finite-thickness, moderate-speed
field benchmark, validated first by local scattering tests and independent
references. Reduced-dimensional fields may help explore the strong limit,
but cannot certify 3D stability, isotropic gas pressure or collision yield.

## Finite-tension algebra must be audited

Using Eq. (5.1), P=dV+2sigma/R, and A=3V/R, the unsimplified energy in
Eq. (5.2) implies

    E = sigma A+(dV+3P)V = dV V(4+3Rc/R),  Rc=3sigma/dV.

The displayed finite-R equality in Eq. (5.2) instead has 4+Rc/R (also
visible in HTML Eq. (11)). These cannot all be exact under the stated
definitions. The bulk limit E≈4dV V is unaffected. We will use the explicit
energy/stress budget and label finite-tension corrections, rather than
silently imposing the printed correction as a numerical success target.

## Next gate

Before any new contracting-pocket production run: validate the new compiled
potential/forces/Hessian against the independent reference; test finite-mass
reflection and frequency shifts with fields; publish the chosen ratios,
stress/spectrum initial data and spatial/time/boundary checks. Contraction,
compression and reversal are observed outcomes. Survival, damping,
equilibrium and cosmological scaling remain separate later tests.
