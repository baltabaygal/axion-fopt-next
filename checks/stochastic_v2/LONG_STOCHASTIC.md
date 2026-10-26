# Longer stochastic evolution: prospective contract

Requested 2026-10-20: evolve the existing stochastic case until the box enters
the true phase. Preserve potential, coupling, backreaction, expansion, periodic
L32, physical seed R8,w2, core tolerance .05, rate Gamma=H0^4 exp(.4*t),H0=.025.
No clipping, extra damping, forced conversion or adjusted nucleation rate.

Full conversion is every sampled lattice cell S>=.5, sustained >=2 elapsed
conformal time units (global histories .04). Separately report all |S-1|<.1,
S range/dispersion and any recrossing. Threshold conversion is not a settled
vacuum. A finite-time equality S=1 is not the completion criterion.

Preserve the original 0..12 candidate file exactly as the prefix. Extend the
Poisson process on disjoint intervals [12,16],[16,18],[18,20],[20,22], then
if necessary [22,24]. Each interval uses the existing sample_candidates,
the unchanged q=Gamma*a^4 sampled with node spacing <=12/1024, and independent
PCG64 seed SeedSequence([0,20261020,index]) with zero-based interval index.
This is a conditional extension of the original realization: rerunning the
old sampler over the whole window would change its past. Save rate tables,
seeds, quadrature errors (<.001), exact candidates and prefix hashes. No extra
false-vacuum factor or geometrical rejection. Freeze the resulting file for
all field controls. No candidates are dropped merely because earlier ones
were rejected. Bound the extension at u24, not an unbounded rate/volume scan.

### Tractability amendment, after assessing u20, before acceleration runs

At u20 only70.5% is converted, with2306 attempts and4 accepts. A certified
rejection shortcut is authorized implementation work to continue this same
physical case. Retain the ordinary solver/build and u20 result. Duplicate the
owned driver for this extension, raise its file budget to1million and extend
the same disjoint-interval scheme to [24,26],[26,28], if needed. Upper bound
u28. Retain every candidate in the input catalogue and log; no rate clipping.

For a base/output interval du<=.02 starting at u, let b=max|S|+1, p=max|P_S|,
a0=a(u),am=a(u+du). With |cos|<=1 the full S-force magnitude is bounded by

    F=lambda/4*(2b+6b²+4b³)+DeltaV*(12b²+12b³)+4f²M²b
    D=am²*12b/dx²+am⁴*F
    B=du*(p+du*D)/a0².

If B<1, a bootstrap bounds every cell's displacement by B through the
interval (also conservatively bounds PV2 full/half kicks and its substeps).
At a candidate, test27 actual grid cells at offsets0,+/-R/(2am) from its
centre, rounded periodically. If any cell has am*distance<=R and
|S_cell|-B>tolerance+1e-10, it remains a rejection witness throughout the
interval. Skip alignment to this candidate, log certified_rejection=1 and
the witness lower bound; exact core maximum/count are NaN/0 because they
were not measured at its time. Otherwise use the ordinary exact-event path.
Recompute after any non-skipped event, integration step or insertion. Enable
only after u12 in physical runs, preserving the old pilot prefix exactly.
Thus the only numerical change is avoiding short steps at certified rejections;
no admissible candidate is deliberately omitted. Energy ledgers remain zero
at rejected events. No simulated energy state is attached to a skipped event.

Validate the new path on a homogeneous converted N4 state and rejected events
against independent DOP853 with du .02,.01,.005: finest error<2e-5, refinement
ratios>3.5, all events rejected and some certified. Rerun the original insertion
tests against the new driver (shortcut inactive there). Compare accelerated
u20 N96 against the saved ordinary run: identical accepted events, endpoint
max component difference<.002 and axion energy relative difference<.001;
retain prefix and all prior accounting gates. Then proceed to longer runs.

The initial homogeneous certificate diagnostic used centre(2,3,4) on a grid
of spacing2: its physical R1.5 core loses every cell as a increases, so most
late candidates followed the ordinary path and refinement ratios3.20/3.37
failed the declared3.5 gate despite tiny errors. Preserve rejection-validation-01.
For the intended continuously exercised certificate control, centre the core
on an actual cell at(3,3,3), retaining the same physics, events and gates.
Report the corrected diagnostic separately; do not overwrite the failed test.

The first accelerated u20 comparison passed acceptance, prefix and energy
(relative3.77e-6), but failed the prospective .002 maximum component gate:
max delta S=.001303,max delta P_S=.013647. Preserve the failed acceleration
audit. Avoiding numerous short rejected-event steps increases the integration
step, so refine the late-time cap instead of loosening the gate. A new build
uses rejection_dt_factor=.25: maximum late du=.005 for the main/base .02 run
and .0025 for its base .01 temporal control, while preserving the original
base timesteps and exact pilot prefix through u12. Repeat both reference
tests and the u20 comparison against the same ordinary reference before
extending farther. The candidate model and numerical gates remain unchanged.

The refined accelerated run differed from the original ordinary run by
.002169 in momentum, slightly failing .002 again, while delta S<.000198
and relative energy difference4.09e-6. Extrapolation of the two accelerated
timesteps shows the original ordinary reference itself has a truncation-error
floor at this tolerance. Preserve this second failed audit. Add a switch to
disable certificates while retaining the refined .005 late timestep, and
run a MATCHED ordinary fine reference with exactly the same early trajectory,
rate and late maximum timestep. Apply the same .002/.001 comparison gates
against that new reference before continuing. This refines the reference;
it does not retroactively make the two coarse-reference audits pass.

### Further duration amendment after u28 assessment

The matched ordinary fine comparison passes: max component .0006335,
relative E_A3.56e-8, identical four accepts and exact original prefix.
The u28 run passes numerical accounting but has no all-cell interval:
98.2924% S>=.5, S_min=-.528, S_max=1.685. Preserve this unmet completion result.
Continue one bounded stage to u40 with the same physical model. Generating
the exponentially increasing candidate process through40 would require about
6.5billion attempted births. Do not clip its rate or stop it by assumption.

Instead certify EMPTY nucleation eligibility throughout EACH step after u28.
Use the same displacement bound B above and divide the lattice into blocks
with b=floor(R/(a_max*sqrt(3)*dx)+.5-epsilon). A candidate centre's nearest
cell is in one block; any witness cell in that block is within
sqrt(3)*(b-.5)*dx of the centre, including its half-cell offset. For every
block, require at least one sampled cell with |S|-B>.05+1e-10. This proves
every possible continuous candidate centre has a forbidden core witness for
the whole step. Periodic distances only reduce this conservative bound.
Log every proof interval, B, b and the minimum witness lower bound. Abort if
any block fails; never silently evolve past an uncertified interval without
explicit candidates. There is then exactly zero accepted intensity, regardless
of Gamma. Omit materializing these provably rejected Poisson marks, not their
physical rate. No random count or positions are claimed for the omitted marks;
report their expected count and the explicit prefix catalogue separately.

The original542713 candidates through28 remain fixed. Independent Python
checks at u24 and28 already give minimum tile witness bounds .827 and .884
for du=.005, well above .05. Add a homogeneous independent ODE control for
the global path at three timesteps, same error/refinement gates. Assert proof
intervals cover28..40 with no gaps, positive duration and every witness>.05.
Apply the original energy, phase and refinement gates at this longer endpoint.
This is analytical rejection on a prescribed numerical trajectory with a
conservative displacement envelope, not a switch-off of nucleation.

### Final bounded duration stage after u40 assessment

u40 has99.5748% S>=.5, S_min=-.194, no all-cell interval. Its2400 exclusion
proofs cover28..40 continuously, minimum witness .5919; numerical gates pass.
Retain the unmet completion result. Extend once to a ceiling u80 (a<=3),
with early stop after THREE sampled conformal units of every-cell S>=.5.
This gives the two controls one extra unit beyond the original >=2 criterion.
Read the actual final snapshot time, not the ceiling, for all diagnostics.
Run controls to that fixed actual endpoint. Keep the same late dt .005/.0025,
all physical inputs, and exclusion certificates. Per-run timeout1800 seconds.
No further duration scan is authorized by this bounded stage. If it remains
incomplete, report persistent remnants and numerical limitations, not forced
conversion or an established stability mechanism. Add a small independent
ODE check of the stop-time/snapshot semantics before using the new stop.

Resource amendment before the two final controls: the independent spatial and
temporal runs may execute concurrently, each with one CPU thread and separate
outputs, at most two native field processes. This changes only scheduling;
there are no comparative performance claims from their wall-clock timings.
The main trajectory still runs first and determines their fixed endpoint.

First assess main N96,maximum du=.02 through u20. If the sustained criterion
is not met, assess u22 and, if needed and tractable, u24. Rerun from the
unchanged initial fields to avoid an unvalidated restart/time-reset. Output
times and candidate prefix are unchanged, permitting overlap checks. The
current backend budget is10000 events; extending it beyond u22 requires an
explicit resource-guard-only change, fresh build and reference validation.

Once the first sustained interval is located, keep that endpoint and frozen
candidate file for N128 at du=.02 (finer spatial control) and N96,du=.01
(temporal control). Sequential one-thread native runs, timeout900 seconds
per run, raw outputs in fresh outputs/stochastic_v2 folders. Histories .04,
slices .2; full initial/final and accepted insertion snapshots. Record the
same acceptance decisions or report any changed catalogues as a limitation.

Numerical gates: finite fields, exact initial transfer, prefix-overlap histories
within1e-10 absolute, independent endpoint energy1e-11 relative to the fixed
max(initial/final total Ebar,sum absolute insertion work) scale, total/axion
expansion+insertion balances<.005 in that normalization. Verify accepted seed
maps (<2e-14), energy jumps (<1e-10 relative), unchanged A/P_A at insertion.
Verify endpoint phase extrema/fractions/moments independently (fraction exact,
moments within1e-12). Where acceptance agrees, endpoint Ebar_A change<10%
spatially and<2% temporally; all controls must show >=2 all-cell interval for
a validated sustained-conversion claim, with first interval start shift<=.4
spatially and<=.08 temporally. If any gate fails, retain and explain it.

Report seed resolution and the later physical spacing a*dx. Colliding,
accelerating walls need their own local resolution study; two grids and phase
convergence are not a continuum wall/spectrum claim. Do not use a spherical
isolated-wall width after interaction. No ensemble, abundance, fully thermal
nucleation, Friedmann/reheating or Marsh exponential-coupling claim.
