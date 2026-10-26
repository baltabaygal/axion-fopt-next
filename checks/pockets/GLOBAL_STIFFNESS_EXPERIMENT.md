# Globally stiffer wall with moving-volume accounting

Declared before runs, 2026-10-17. Previous turn completed the structured-tail
test: b=12 worsens the axion-sector survival screen. Raising vacuum curvature
alone does not address the large forcing near the barrier core. This test
raises the whole barrier, explicitly changing tension and bound spectrum.
It cannot isolate one leakage mechanism or reproduce the paper's potential.

Use static coupled fields, p=4, b=0, lambda=32 (tanh width w=.5),
DeltaV=.01, R0=20, m_true=2, seed 0, .1<=omega<=.6. Degenerate tension
sigma=2/3, false/true dilaton masses 4 and sqrt(16.12). Thin-wall mean
pressure balance gives Ea=3(DeltaV+2sigma/R0)*4piR0^3/3=7707.373976806957.
This is an initial estimate, not measured equilibrium. Mean required pressure
divided by degenerate barrier height is about .153 (previous w=1: .347).
Local coherent stresses can still be large. Modes are Gaussian classical
equipartition, not a demonstrated thermal population.

Sample the same radial eigenfunctions (h=.01,rmax=31,rcut=30), seed and
amplitude on N=320 and 384, L=64. Fix amplitude only on N=320.
The initial 10--90 wall width is 1.09861: 5.49 and 6.59 cells.
Maximum-frequency exterior axion penetration is .52414: 2.62 and 3.14
cells. R0/w=40. These are finite hierarchies, not the bulk cosmological
limit. Use dt=.0125, t=80, histories .25, slices 4, full states 40.
dt*m_true=.025, dt*m_phi,TV=.0502; the driver also guards the actual
Hessian/lattice frequency bound. Separate timestep refinement remains due.

Enable fixed radius 26 and the independently validated moving monitor
Rin=26,Rout=30. Require the false-vacuum body to remain inside Rin at all
saved histories; all moving lab/advective terms remain in raw output.
This window is not an absorbing boundary; periodic returns are possible.
Moving storage/net current is not irreversible escape or invariant mass.

Run sequentially with four OpenMP threads. Use the streamed I/O build
validated in `outputs/pockets/validation-tracking-02`; archive source/build
hashes before launch. A large-grid memory failure invalidates that run.

Gates retain prior screens: exact transfer, completion, global drift<5e-4;
fixed and moving total/sector continuity residuals and twice-coarser cadence
differences <=1e-3 of initial respective window energy; moving enclosure;
all slice and final 3D wall widths measured unambiguously and >=4 cells.
Paired radius histories must differ <=.25, far-axion fractions <=.005,
final sector depletion <=.01. Candidate additionally requires R>=16 and
axion-sector depletion<=5%; report total moving retention independently.
Do not retroactively change older failed classifications using this monitor.

A passing pilot only permits time, box and perturbation checks followed by
t>=400 evolution. It establishes no long-lived relic, mass-radius law,
expansion survival, formation probability or physical mass/radius scale.
