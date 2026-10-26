# Instructions for the axion relic pocket workspace

Read README.md here, ../AGENTS.md and ../../docs/physics.md before changing
this workspace. It owns the axion relic pocket (ARP) experiments and is
separate from the homogeneous, single-bubble, stochastic and PQ-string
workspaces. The reference paper is 2609.24714.pdf in this directory. It is
read-only; cite its equation numbers when implementing a formula.

- Separate what the reference paper derives (trapped-mode photon conversion
  and phenomenology) from what its Ref. [5] (arXiv:2407.01676) derives
  (formation, R_eq, epsilon, alpha). The 2024 formation paper contains a
  spherical reflective-particle simulation (Sec. 5, Fig. 4), not a coupled
  field simulation; it does not establish the damping timescale. The 2026
  phenomenology paper does not contain field simulations. New work follows
  the user-approved FIELD_GOAL_REVISED.md and remains fields only.
- The paper's numbers are physical (GeV, metres) and our lattice units are
  dimensionless. Never infer a physical pocket mass, radius or transition
  temperature from a toy run. State the dimensionless hierarchy actually
  simulated: pocket radius, wall width, axion wavelengths, m_true and dx.
- Reflection requires omega < m_true. Report the penetration depth
  1/sqrt(m_true^2 - omega^2) and wall width in cells, and treat the Dirichlet
  boundary (3.2) as a limit to approach, not an assumption. Check that dt*m_true
  and the largest lattice frequency are resolved.
- The formation paper also discusses nonthermal axions; thermalisation is
  not imposed on our initial data or claimed without evidence. A classical field
  represents only high-occupation modes. Record the spectrum, amplitude,
  occupation assumption, filter and seed of any random axion initial state.
  Distinguish it from coherent misalignment initial data.
- Prescribed arrival-map walls cannot be stopped by axion pressure. Do not
  claim pocket stabilisation, R_eq, epsilon or alpha from prescribed geometry.
  Those need a dynamical wall (a dilaton field coupled to the axion mass) with
  an explicit, documented potential and coupling.
- Use the actual compiled CosmoLattice adapter for lattice results. Keep the
  existing adapter, its tests and the other workspaces working. New
  model/driver code belongs here; builds go under the ignored root build/.
  Never edit upstream vendor source.
- Compare with independent references: spherical Bessel zeros for trapped
  modes, an independent radial/ODE solver for finite-m_true leakage, and the
  paper's overlap formulas (3.30)-(3.36) for any emitted-power estimate.
- Predeclare numerical gates before each comparison. Refine time and space
  separately with fixed physical parameters. Hold geometry or initial states
  fixed when comparing discretisations.
- Save parameters, initial states, seeds, build identity, histories, spectra,
  references and machine-readable checks. Write raw outputs into fresh
  directories under outputs/pockets/. Compact evidence may live in results/
  here. Never overwrite earlier runs.
- Keep runs modest and report what each check establishes. Physical
  detector rates, dark-matter abundance claims and manuscripts are separate
  tasks.
