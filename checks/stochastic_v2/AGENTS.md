# Instructions for stochastic_v2

Read README.md here, ../AGENTS.md, ../../docs/physics.md and the reference
provenance in references/README.md before doing work in this folder.

## Current authorization

On 2026-10-20 the user requested longer stochastic evolution until the full
box is in the true phase. Read LONG_STOCHASTIC.md for its prospective gates
and recorded amendments. The u20/u24/u28/u40 stages all have no full-cell
interval; retain those unmet completion results. A final bounded stage to
u80 with early stop after three sampled all-cell S>=.5 units is complete,
including spatial/time controls at u80. No
changed friction, clipping, forced seeds or physical parameters are authorized.
The long driver certifies rejected candidates and then globally empty seed
eligibility after28; it aborts if that sufficient proof fails. Do not describe
this as shutting off the physical rate. Read the contract before follow-up.

Read results/long-stochastic-01/README.md and report.json for the final long
checkpoint. All numerical gates pass, but the persistent-conversion goal is
UNMET: no sampled all-cell interval exists in any of the three runs. Main N96
has342 below-half cells at80; N128 has261, and half-dt N96 has346. Axion-energy
differences4.23% spatial/.00309% temporal pass the declared pilot gates.
The remaining physical volume is strongly grid dependent (342 versus110.1),
so do not claim converged pocket counts, volume, stability or an axion-pressure
cause. Main audited replayed full checkpoints show every cell exceeded.500001
at least once; all342 endpoint below-threshold cells are recrossings. That is
not simultaneous or sustained conversion, an exact first-passage time, or a
refined coverage-map result. Preserve the explicit goal_reached=false report,
all failed earlier comparisons and intermediate unmet duration results. No
further duration, resolution or physical-model scan is automatically requested.

The user explicitly requested an active goal to add expansion, validate a
coupled-field single bubble, then implement/validate a small stochastic
nucleation pilot. That implementation and controlled pilot are complete.
Prescribed radiation domination and
phenomenological finite-radius insertions are in scope; retain polynomial
mass coupling. Large ensembles, microphysical tunnelling/thermal models and
switching to the Marsh exponential coupling are separate physical requests.

Read COSMOLOGY_INPUTS.md, EXPANSION_STOCHASTIC_CONTRACT.md and
results/expanding-single-bubble-01/README.md plus
results/expanding-stochastic-01/README.md for the new checkpoint. The expanding
single bubble and stochastic seed-0 pilot passed their declared gates. Keep
expansion and insertion work in energy balances; snapshots store u and P=a²X'.
The frozen 12-candidate file produced three accepted seeds and nine rejections,
identically in the two controls. This is one short realization, not completion
or abundance evidence. No additional runs or production scans are automatically
authorized by this completed checkpoint. The older initial insertion build
failure and aborted empty-default control remain recorded in outputs/build;
successful validation uses insertion build03 and insertion-validation-02.

The user initially requested research only. On 2026-10-19 the user said
"lets do that first experiment", authorizing implementation and controlled
single-bubble runs, including the prerequisite validation. Large stochastic
ensembles and production scans remain outside this first-experiment scope.
Do not repeatedly request approval for already authorized research/setup work.

## Scientific contract

- Distinguish the existing dynamical axion from the proposed transition-driving
  scalar S. A smooth prescribed mass wall is not a dynamical order parameter.
- Mark proposals, assumptions, open questions and measured findings separately.
  README equations and candidate potentials are not selected production physics.
- Specify the full action, kinetic normalization, U(S,T), m_a^2(S,T), units,
  background and bath before future implementation. Retain all S-dependent
  offsets when matching the effective potential to Appendix A.
- Include the reciprocal axion force on S for a coupled model. If omitted,
  declare and justify the spectator approximation. Check the full potential's
  vacua and barrier at the actual axion initial angle.
- Classical evolution does not itself supply metastable-vacuum tunnelling.
  Separate prescribed bubble insertion, bounce-derived rates and thermal noise.
  Do not describe seeded real-time field evolution as spontaneous nucleation.
- Do not impose the old constant-speed geometry identities on dynamical walls.
  Finite birth radius, wall acceleration, collisions and recrossing need new
  contracts. One first-arrival map cannot replace the evolving S configuration.
- A slow terminal wall speed requires a friction/transport model. Distinguish
  phenomenological damping from hydrodynamics, reheating and a self-consistent bath.
- State cosmic/conformal time and physical/comoving lengths. Derive scale-factor
  and field-rescaling factors explicitly. Keep f^2 in the S feedback force.
- Account for inter-sector exchange, expansion, explicit time/temperature
  dependence, damping and bubble-insertion sources separately. Do not require
  closed-system energy conservation with unaccounted external sources.
- Future numerical evidence must use the actual compiled field backend plus
  independent references. Predeclare gates and refine time, wall/field spatial
  scales, boundaries, volume and ensembles separately as appropriate.
- Audit scalar/axion scale separation and Lorentz-contracted wall resolution
  before choosing a grid. A small unresolved wall is not a thin-wall validation.
- Start future implementation with homogeneous and single-wall/bubble controls;
  do not jump directly to full 3D stochastic production or abundance claims.

## Repository boundaries and provenance

- Preserve checks/stochastic and other workspaces as references. Do not edit
  their results, third-party sources or unrelated local work during setup.
- Inspect checks/pockets for reusable coupled-field and accounting methods before
  rebuilding them. Its nonzero false-vacuum mass/exponential coupling and pocket
  physics do not establish this workspace's massless false-vacuum model.
- After implementation is requested, use fresh outputs/stochastic_v2/ run folders
  for generated outputs. Keep compact reviewed evidence here. Save source/build
  provenance, parameters, seeds, insertion events, raw data and failed gates.
- Keep the pinned paper PDF and text intact. Download later versions under new
  versioned filenames and record URL, retrieval date and SHA-256.
- No delegation is requested by these instructions.

## Current checkpoint

The first experiment is complete. The dedicated owned model uses a massless
false-vacuum cosine axion and full reciprocal coupling, with an actual 3D
CosmoLattice PV2 backend. It does not modify pockets or the original adapter.
Read EXPERIMENT.md and results/single-bubble-02/README.md for parameters,
validation, limitations and exact reproduction commands. Retain the failed
initial configuration attempt and successful reference evidence. The pilot
does not establish a continuum limit, terminal wall speed or cosmological yield.

The user subsequently requested longer single-bubble evolution and every input
parameter. That extension is complete through t=18; see LONG_EXPERIMENT.md,
INPUT_PARAMETERS.md and results/single-bubble-long-01/README.md. The current
backend uses compensated host-side diagnostic reductions and supports N<=192.
The original large-box snapshot accounting failed a 1e-11 gate; retain that
report and the passing compensated replay audit. Do not overwrite the original
histories or represent all original gates as having passed. No longer run or
stochastic production is automatically requested by this checkpoint.

The user then requested a smaller box to obtain full conversion. L24 and L20
were simulated, with phase statistics added passively to the owned driver.
Read results/periodic-conversion-01/README.md and the three SMALL_BOX*.md
contracts before follow-up. L20 has a refined, measured transient all-cell
S>=.5 interval near t28..29, but no >=2 sustained interval and no relaxed
uniform S=1. Preserve the L24 dt=.02 energy failures, the original L20
inapplicable spherical-width failure, and the separately passing transient
refinement report. Do not force conversion, clip excursions or add damping
without a new physical request/contract. No unbounded smaller-box scan or
further simulation is automatically requested by this checkpoint.
