# Expansion and stochastic-pilot inputs

For the 2026-10-20 duration extension, use
[LONG_STOCHASTIC.md](LONG_STOCHASTIC.md): it preserves this physical model but
adds conditional Poisson intervals through28, certified empty eligibility
after28, a smaller late timestep, and a bounded duration up to80. The short
pilot table below is historical; it does not describe the longer endpoint.

These are the implemented 2026-10-19 extensions. The earlier static runs keep
their own inputs in INPUT_PARAMETERS.md. No parameter was fitted to the paper.
The prospective tests and thresholds are in EXPANSION_STOCHASTIC_CONTRACT.md.

## Fields and background

Natural units, v=c=hbar=1; canonical S and A=f theta. Length/time units v^-1,
fields v, energy density v^4. lambda=2, DeltaV=.06, M=.5, f=.2, theta_i=1.

    V = lambda*S²*(1-S)²/4 + DeltaV*(1-4*S³+3*S⁴)
        + f²*M²*S²*[1-cos(A/f)]
    a(u)=1+H0*u, H0=.025
    elapsed cosmic t=u+H0*u²/2
    Hphysical=H0/a², Hconformal=H0/a
    P=a² X_prime=a³ X_dot
    X_prime=P/a²; P_prime=a² laplacian(X)-a⁴ V_X

Both fields evolve and the axion exerts its reciprocal force on S. This keeps
the polynomial mass coupling, with zero axion mass at S=0. It does not adopt
the Marsh/pockets exponential mass. Wall velocity is an output. Friction,
thermal noise, temperature dependence, tensor modes and a dynamical bath are
absent. The prescribed RD background is not a solved Friedmann equation:
subdominance relative to a radiation density requires a future physical scale
map, and is not established by these dimensionless runs.

## Runs

| Setting | Expanding single bubble | Stochastic pilot |
|---|---|---|
| Comoving periodic box L | 24 | 32 |
| Main grid | 96³ | 96³ |
| Spatial control | 64³, same dt | 64³, same dt and candidate file |
| Main conformal dt | .02 | .02 maximum; short event/output steps |
| Temporal control | 96³, dt=.01 | 96³, maximum dt=.01, same candidate file |
| Conformal interval u | 0..16 | 0..12 |
| Cosmic elapsed endpoint | 19.2 | 13.8 |
| Scale factor | 1..1.4 | 1..1.3 |
| Initial S | tanh bubble R=8,w=2 at box centre | exactly zero |
| Initial A | uniform .2 | uniform .2 |
| Initial P_S,P_A | zero | zero |
| History cadence in u | .04 | .04 plus pre/post insertion rows |
| Central slice cadence in u | .2 | .2 |
| Full states | initial/final | initial/final and each accepted pre/post seed |

Cell centres (i+.5)L/N, three dimensions, periodic all axes, double precision,
nearest-neighbor Laplacian, compatible forward-difference gradient energy.
Actual CosmoLattice Position Verlet PV2; one CPU process/thread, C++20 -O2,
no fast-math. Compensated host energy sums. Framework alpha=fStar=omegaStar=1,
omegaEoS=1/3, a0=1,t0=0, expansion=true,fixedBackground=true, withGWs=false,
PS_type=1; kIR=2*pi/L follows from lSide. No field rescaling is applied.
OMP_NUM_THREADS=1, OMP_DYNAMIC=FALSE, OMP_PROC_BIND=false,
OPENBLAS_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1.

## Stochastic process and seed map

Physical rate Gamma=H0^4 exp(beta*t_elapsed), beta=.4 (beta/H0=16),
Gamma0=3.90625e-7. Candidate intensity per comoving volume/conformal time is
q=Gamma*a^4. The existing core.sample_candidates samples the linear interpolant
of a saved 1025-node rate table on u0..12, uniform comoving positions, seed 0,
NumPy default_rng/PCG64. Candidate budget 10000. Its wall_speed=1 is unused
catalogue metadata, not a prescribed scalar-wall speed. The 12-candidate file
is frozen and copied byte-for-byte for controls. No extra false-vacuum factor
and no old arrival-map/geometrical exclusion are used.

At an exact candidate time, all cell centres inside physical radius Rseed=8
must have |S|<=.05. At least one core cell is required. Physical r uses a times
the nearest-image comoving distance. The profile b is the tanh bubble with
w=2, tapered with 1-10*z³+15*z⁴-6*z⁵ from r=12 to16, then exactly zero.
The support must fit inside half the physical box; the implementation checks
this conservatively at a=1. For accepted candidates:

    S_new=S_old+(1-S_old)*b; P_S,new=(1-b)*P_S,old
    A_new=A_old; P_A,new=P_A,old.

The seed core is approximately at rest in the comoving bath; its tanh b is
not identically one. This finite-core exclusion and finite seed radius are
phenomenological modeling choices. The rate and profiles are not computed
from a bounce or a bath. The initial tanh 10--90% physical width is
2*w*atanh(.8)=4.394449; seed-width cell counts are reported at each accepted
birth. Those counts do not establish resolution of later accelerating or
colliding walls. Multiple simultaneous candidates are processed in file order.

## Energy convention

Histories report physical densities integrated over COMOVING volume, Ebar;
physical total box energy is a³ Ebar. K uses P²/(2a⁶), G uses (grad X)²/(2a²).
Ebar_A begins at exactly zero in the stochastic run because S=P_A=grad A=0.
Accepted insertions change the S-dependent axion potential immediately even
though A and P_A are unchanged. They are external work, recorded as I_total,
I_S and I_A, with no assumption that every jump is positive.

    Ebar_total'= -Hconformal*(6K_total+2G_total) + insertion impulses
    Ebar_A'= Q_A-Hconformal*(6K_A+2G_A) + axion insertion impulses
    Q_A=integral 2*f²*M²*S*(P_S/a²)*[1-cos(A/f)] d³x

The balances integrate between pre/post event rows at identical u and subtract
the exact jump ledger. Use physical scale factors when interpreting snapshots:
their legacy time slot stores elapsed conformal u, and their momenta are P.
Each run saves input.in, cosmology.json, build hashes, histories and snapshots.
