# Field-only relic pocket programme

**Historical programme.** The user accepted a revised scientific goal on
2026-10-18: [FIELD_GOAL_REVISED.md](FIELD_GOAL_REVISED.md). It now governs
new work. The original sequence and evidence below are retained intact.

Agreed 2026-10-17. Use coupled classical axion and dilaton fields only.
No particle or N-body implementation is part of this programme.

## Success sequence

1. Establish a pressure-supported static-background pocket with resolved
   walls, bounded size under small perturbations, and a converged small loss
   rate over at least ten initial diameter-crossing times (t >= 20 R0).
   Survival without a controlled loss rate is not sufficient.
2. Across several radii and energies, test pressure balance and mass-radius
   scaling, with finite tension, finite mass, non-sphericity and radiation
   corrections. Physical mass/radius conversions require an explicit choice
   of physical scales; formation has to supply the initial size distribution.
3. Extend the validated coupled field model to prescribed radiation-dominated
   expansion. Vary H_* R0 and measure physical radius, retained energy and
   lifetime in Hubble times. Use continuity, not constant total energy, as
   the expansion conservation gate.
4. Study field pockets produced by bubble collisions, including initial
   velocities, sizes, formation efficiency and population dilution.

The static pocket and expanding survival studies come first. Do not infer a
formation efficiency or cosmological abundance from an imposed initial ball.

## Current evidence and limitations

Prototype-01 and leak-02 use m_a^2(phi)=m_TV^2 phi^2. Their gas is a seeded,
high-occupation classical field with a finite mode band, not a Bose thermal
ensemble. The initial wall is a tanh profile, not a self-consistent wall-gas
solution. Existing leak runs demonstrate delayed collapse, not stability.

`audit_leak.py` reads a snapshot of an existing suite into a fresh directory.
It never changes the original run. It distinguishes unfinished histories,
missing samples and a vanished pocket. Far-field stored energy is not
cumulative loss; periodic radiation can return. For L=64, R0=20, the initial
surface-to-image-surface light travel distance is only 24. This is an
isolation warning, not a measured return time. A larger box or a separately
validated treatment of outgoing radiation is required for lifetime claims.

## Next pilot, declared before launch

Keep the compiled phi^2 coupling unchanged, and test whether removing the
large initial pressure deficit reduces the violent transient. This is a
diagnostic pilot, not the ten-crossing-time success experiment.

Use the existing bound-mode field construction: R0=20, w=2, lambda=2,
DeltaV=0.01, m_TV=2, L=64; 0.15 <= omega <= 0.6, radial h=0.01,
r_max=31, r_cut=30, equipartition Gaussian coefficients, seed=0.
Assume only high-occupation classical modes. No thermalisation is imposed.
Set the initial lattice axion energy to

    E_a = 3 (DeltaV + 2 sigma/R0) (4 pi/3) R0^3,
    sigma = sqrt(lambda/2)/6.

This is a thin-wall isotropic pressure estimate, following the equilibrium
condition in arXiv:2407.01676 Eq. (10); it is not an assertion that the
inhomogeneous field stress is locally balanced. Reference paper 2609.24714
Eq. (2.3) uses the associated equilibrium size approximation.

N=128, dx=0.5, dt=0.0125, duration=80 (two initial diameter crossings).
Record histories every 0.5, slices every 4, full initial/final states.
The exterior penetration depth at omega=0.6 is 0.524 (1.05 cells);
the tanh wall's initial 10-90 thickness is 4.39 (8.79 cells).
Finite-mass evanescent structure therefore still needs refinement.

Pilot numerical gates: exact state transfer; energy drift <5e-4;
10-90 wall width >=4 cells in saved central-row slices while R>=8;
no guard stop before t=80. Report all failures. A candidate is worth
refining if R remains >=16 to t=80 and the stored far axion energy at t=80
is <5% of initial E_a. These are screening thresholds only: the periodic
box is not isolated over this duration. Do not call this relic survival.

Next decisions: compare with existing low-pressure initial data, resolve
the tail/wall and box effects, then consider a documented alternative mass
function whose false-vacuum derivatives do not strongly shift the dilaton
mass. Do not change the potential silently or fit a lifetime by tuning it.

## Progress after the first continuation

The balanced initial state passes early history comparisons for timestep,
space and box size; it still has thin/ambiguous wall regions and transfers
about 24% of axion-sector energy into the dilaton by t=80. Relaxing the
mean wall at fixed classical actions does not remove this depletion.
See [results/balanced-refinement-01.md](results/balanced-refinement-01.md).

A separately documented phi^4 mass-function comparison reduces depletion
to 13% but fails its factor-two improvement threshold and slice diagnostic.
The default phi^2 evolution is reproduced. See
[results/coupling-pilot-01.md](results/coupling-pilot-01.md).
The extended compiled model passes 40 validation checks.

No long-lived, converged field pocket has yet been established. The next
controlled question is whether separating the axion frequencies from the
dilaton radiation threshold reduces the transfer. Keep all original
success requirements and the expansion/formation stages intact.

## Frequency hierarchy and energy-current diagnostics

The three-band frequency pilot also fails the <=5% depletion candidate
criterion: depletion is 13.2%, 9.4% and 11.5% for soft/intermediate/broad
states. Soft states have fewer modes and greater asphericity; changing the
initial pair-frequency ceiling alone does not isolate the radiation cause.
See [results/frequency-pilot-01.md](results/frequency-pilot-01.md).

An optional fixed-sphere symmetric-bond energy monitor now measures
outgoing/incoming net currents and the coupling exchange between fields.
It changes no force or boundary. The compiled monitor matches an independent
Python implementation and its continuity residual converges at second
order for both mass functions; 51 compiled validation checks pass.
Its physical-run output cadence must still be checked.

The completed suite `outputs/pockets/stiff-wall-01`, declared in
[STIFF_WALL_EXPERIMENT.md](STIFF_WALL_EXPERIMENT.md). It doubles the
dilaton mass (and tension), resolves the narrower wall on N=192 and 256,
holds one continuum gas amplitude fixed between grids, and records the new
energy-current monitor. Its histories agree between grids, but depletion
and wall diagnostics fail the original candidate gates. See
[results/stiff-wall-01.md](results/stiff-wall-01.md).

## Latest evidence and remaining work

The fine pocket drifts 6.07 units by t=80 and partly leaves the fixed monitor.
The fixed-sphere export cannot be equated with radiation loss. A periodic
centroid snapshot audit encloses the body and retains 97.4% of initial stored
energy at t=80; it is not a moving-volume continuity or ten-crossing-time
result. Quantifying body motion and radiation separately remains necessary.

Compensated diagnostic sums correct a tight initial-energy accounting
roundoff failure without changing any field evolution or relaxing the old
threshold. The old failed reports are retained. Optional OpenMP builds pass
the same 51 pocket checks and the main-adapter checks; matching fields and
static checkpoint continuation are bitwise identical in controlled tests.
These are supporting numerical capabilities, not goal completion.

The completed `outputs/pockets/empty-stiff-control-01` checks the role
of gas for the stiff potential in a trusted early contraction window.
All early gates pass; the empty reference reaches R=15 by t=9.75 while
the matched gas pocket remains at R=19.469. This is pressure-support
evidence, not a long lifetime. See
[results/empty-stiff-control-01.md](results/empty-stiff-control-01.md).
The independent structured-wall potential/kink reference and compiled
implementation pass their gates, including independent ODE and continuity
refinement. The b=0 implementation reproduces the prior fields exactly.
The physical paired-grid suite is complete. Both b=12 cases fail the original
survival screen, with 19.3% axion-sector depletion and failed wall diagnostics.
Paired histories agree, and a moving-centroid snapshot audit retains 94.3%
of total stored energy at t=80. See
[results/structured-wall-01.md](results/structured-wall-01.md).

A moving weighted-volume monitor now passes independent differential
continuity and periodic-translation tests and compiled reference checks
(`outputs/pockets/validation-tracking-02`, including streamed I/O changes). It subtracts boundary advection
from lab-frame currents without changing field evolution. These checks do
not yet establish physical-run retention. Streamed binary I/O and prepared
launches pass format/corruption tests; the enlarged grid cap is numerical
capacity only, with no large-grid physics result yet. All original
survival, mass-radius, expansion and formation requirements remain active.

The next declared pilot is [GLOBAL_STIFFNESS_EXPERIMENT.md](GLOBAL_STIFFNESS_EXPERIMENT.md):
lambda=32, w=.5 at N=320/384, with fixed continuum amplitude and both fixed
and moving energy monitors. `outputs/pockets/global-stiffness-01` holds its
source/build contract and raw results as they are produced. The runner's
moving-current path passes a compiled smoke check in
`outputs/pockets/tracking-runner-smoke-01`; repository tests pass (72,
11 optional skips). This is numerical readiness, not survival evidence.
The lower-memory full-state wall audit reproduces the saved N=256 result
exactly (`outputs/pockets/wall-memory-equivalence-01`), including every
ambiguous contour point and its failed classification. The running global
stiffness suite keeps its original loaded audit; new analysis uses fresh
outputs. Repository tests after this change pass (74, 11 optional skips).
The N=320 globally stiff run has entered field evolution. Its exact initial
transfer, independent energy accounting, frequency guard and initial slice
wall check pass; see [the initial audit](results/global-stiffness-initial-01.md).
The fine case and physical survival checks remain pending.
At t=40 the coarse globally stiff pocket retains 99.925% of moving-window
energy and its radius has rebounded after contraction. Independent snapshot
energy accounting passes, but the full-3D wall has 261 of 192736 edges below
four cells (minimum 3.63): the coarse midpoint resolution check fails.
See [the midpoint evidence](results/global-stiffness-initial-01.md).
The finer run remains required; no lifetime or scaling claim follows.
The coarse globally stiff case has completed t=80 and fails its original
sector-depletion and final wall-width screens; the finer state is preparing.
Moving net export rises strongly over the four successive 20-unit intervals.
Additional [moving-grid comparison gates](MOVING_GRID_COMPARISON.md) were
declared before finer-state/evolution results were available. They compare
normalized storage, integrated net current/exchange and periodic centroid
motion at shared times, while retaining all original failed classifications.

The globally stiff N=384 case has now completed. Both paired-history and
moving-current grid comparisons pass, but both candidate screens fail:
fine axion depletion is 6.914%, minimum final wall width 2.655 cells.
Fine moving-window energy retention is 98.780% at t=80. See
[the completed pair](results/global-stiffness-pair-01.md). This supersedes
earlier pending-run status; the full goal remains unachieved.
The thinnest patches correlate with rapid level-set motion, but some
inferred speeds exceed one and cannot certify a material-wall interpretation.
Further uniform refinement exceeds current memory capacity on the verified
16-GB host. The first [memory reduction](LOW_MEMORY_CONTRACT.md) replaces
cubic tracking windows with plane buffers and passes compiled validation
plus bitwise moving-pocket comparisons. No grid-cap increase or new
large physics run has been made; redundant input/state copies remain next.
The diagnostic full-state copy is now replaced with const field views and
chunked snapshot writes. Compiled reference checks and bitwise full/slice/
frame/history comparisons pass (`validation-stateview-01`,
`stateview-agreement-01`). The modest direct-read timing is slower, not a
claimed speedup. Initial-state duplication and a large-grid memory audit
remain before further refinement; the cap is still 384.
