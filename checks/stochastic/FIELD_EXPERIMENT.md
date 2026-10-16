# First stochastic axion evolution (before execution)

Use the indexed pinned-paper tracer driver, 100000 tracers, seed 0,
H=beta=1, L=3, v=0.5, nucleation tau=0..2.5, outer tracer dt=0.05.
JIT warm-up excluded from catalogue timing. Compare the complete catalogue
exactly to the original all-pairs driver. Freeze it for every field run.

Prescribed radiation expansion: a=tau, cosmic t=tau^2/2. Field initial
conditions are imposed at tau=1 (t=0.5,a=1): homogeneous theta=1 and
theta'=0. Earlier bubbles remain in the catalogue; this is not an evolution
of the axion from the Big Bang. End tau=2.5,t=3.125, solver u=tau-1.
Cosine V=m^2(1-cos theta), mf=0, mt=f_axion=omega_star=1.
Use the existing smooth first-arrival mass profile
m^2=0.5[1+tanh((tau-arrival)*v/0.15)] (comoving width 0.15),
including its precursor tails. Geometric phase is the sharp arrival condition,
not a threshold on smooth m^2 or an active-tracer completion flag.
Bubble growth/nucleation is prescribed without axion/plasma backreaction.

Run actual CosmoLattice PV2 at N=32,64 with fixed du=0.005 (300 steps),
and N=64,du=0.0025 (600 steps). Save initial, midpoint u=0.75 and final
binary states, histories and native phase timings. Project exact periodic
arrival maps separately and report their cost, not as tracer engulfment.
The wall spans 1.6 and 3.2 cells, respectively; this is a first resolution
check, not proof of convergence in the continuum.

Independent small-grid check N=8 with du=0.005,0.0025: integrate the same
semi-discrete cosine PDE with DOP853, theta'=pi/a^2,
pi'=a^2*laplacian(theta)-a^4*m^2*sin(theta), periodic nearest-neighbor
stencil. Reference rtol=2e-10,atol=2e-12,max_step=0.01.

Acceptance: exact baseline/indexed catalogue; exact arrival transfer;
all snapshot mass cells agree within 3e-14, a/time mapping within 1e-12;
finite fields/history; snapshot moments/potential match history within 1e-12;
N8 reference max field error <2e-4, momentum error <1e-3, refinement
reduces field error by >3; N64 timestep final-field RMS change <0.001;
N32 vs block-averaged N64 final-field RMS change <0.02 (theta_initial=1).
Record failures without changing gates. Geometry F vs analytic is a
single-realization diagnostic, not an ensemble agreement gate. Do not claim
energy conservation for expansion/prescribed changing mass or final abundance.
