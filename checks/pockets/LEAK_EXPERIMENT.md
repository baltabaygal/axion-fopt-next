# Prototype-02: pocket lifetime and axion leak rate

Recorded before any run. Prototype-01 ([results](results/README.md)) showed
that trapped axions keep the false vacuum open far beyond the no-gas collapse
time, but the pocket leaks and shrinks. A classical field has no Boltzmann
suppression of modes above m_TV, so permanent stability is not a meaningful
target. This experiment instead measures lifetime and leak rate, and asks two
questions:
- Is the leak set by the physics of this toy or by the lattice cutoff?
- How does it depend on m_TV/omega?

As before, the run starts from a ready-made pocket. It does not test formation.

## Model and fixed parameters

Model, units and lattice conventions are as in [EXPERIMENT.md](EXPERIMENT.md):
- L=64;
- wall width w=2 (lambda=2, sigma=1/6);
- dV=0.01;
- R0=20 at the box centre, pi_phi=0;
- static background, PV2.

The changes from prototype-01 are:
- **Timestep** dt=0.025, so the frequency guard bound is about 0.27 at N=128.
- **Gas initial state.** The gas is made of bound pocket modes R_nl(r) Y_lm.
  These are solved for the static initial wall phi(r) from
  -u''+[l(l+1)/r^2+m_TV^2 phi^2]u=omega^2 u by finite differences, with
  h=0.01, u(0)=u(31)=0 and 0.15<=omega<=0.6. Amplitudes are equipartition
  (field variance 1/omega^2, momentum variance 1), seed 0, drawn in
  (l, n, m) order, and evaluated for r<30. The modes fill the pocket up to its
  turning surfaces, with no window gap. Their Dirichlet limit is
  arXiv:2609.24714 Eq. (3.4).
- **Gas energy.** Lattice axion energy, including coupling, equal to
  0.0464*(4 pi/3)*R0^3 = 1554, the same as prototype-01.
- **Far-field measure.** Axion energy beyond R_eff+6 from the false-vacuum
  centroid (periodic nearest image), recorded every 0.5 time units.

## Runs

| Name | N (dx) | dt | m_TV | m_TV/omega_max | Duration |
|---|---|---|---|---|---|
| control | 128 (0.5) | 0.025 | 2, no gas | - | 30 |
| m2_N128 | 128 (0.5) | 0.025 | 2 | 3.33 | 400 |
| m2_N128_dt | 128 (0.5) | 0.0125 | 2 | 3.33 | 100 |
| m2_N96 | 96 (0.667) | 0.025 | 2 | 3.33 | 400 |
| m2_N192 | 192 (0.333) | 0.025 | 2 | 3.33 | 400 |
| m1.5_N128 | 128 | 0.025 | 1.5 | 2.5 | 400 |
| m3_N128 | 128 | 0.025 | 3 | 5 | 400 |
| m4_N128 | 128 | 0.025 | 4 | 6.67 | 400 |

The m_TV=2 runs share one continuum initial state. Its amplitude is fixed
on N=128 and sampled on every grid. The m_TV scan solves its own pocket modes
for each mass, which changes the number of bound states (93 at m_TV<=2, 69 at
3 and 4). This confound is recorded, not removed.

Interpretation note (added after launch, before any result was examined;
no criterion changed): with m_a^2 = m_TV^2 phi^2, raising m_TV changes more
than the trapping margin. Total radiation pressure per reflection is fixed by
momentum transfer, but three other things change with m_TV:
- the interior dilaton mass shift m_TV^2 <a^2>, 4x at m_TV=4 compared with 2;
- the turning surface and force layer, which sit at phi = omega/m_TV;
- the bound-state count.
The m_TV classification therefore describes this coupling family, not the
trapping margin alone. A planned follow-up uses m_a^2 = m_TV^2 phi^4, whose
interior coupling vanishes at phi=0 for every m_TV.

## Measures

- f_far(t) and f_tv(t): axion energy in the far field and in deep true
  vacuum (phi>0.99), each divided by the initial axion energy.
- Leak rate Gamma: least-squares slope of f_far(t) over 25<=t<=125. If the
  pocket dies earlier, the window ends at the lifetime; it needs at least 25
  time units.
- Lifetime t_life: the first t with R<8 (0.4 R0). It is a lower bound if the
  pocket survives the run.
- Control collapse time t_c: the first t with false-vacuum volume below 1% of
  its initial value.
- Late dR/dt: slope of R over the last 100 time units the pocket is alive.

## Numerical gates

1. validate.py passes. This includes a Bessel-limit check of the radial mode
   solver (step wall, mass 400, h=0.01): relative error below 2e-3 for l,n<=2.
   It also includes a far-field consistency check between the compiled driver
   and the Python measure.
2. Initial confinement: f_far(0) < 1e-6 in every gas run.
3. Exact initial-state transfer.
4. Energy drift below 5e-4 while R>=8. With a 10% leak, this keeps the energy
   error under 1% of the leaked energy.
5. No solver guard stop while R>=8. A stop after the pocket dies, or in the
   control after collapse, is an expected end.
6. Timestep, t<=100: max |R(0.025)-R(0.0125)| < 0.05 and
   |f_far difference at t=100| < 0.005.
7. The N=192 wall stays at least 4 cells thick (10-90%, central slice row)
   while R>=8.

## Classifications (reported either way)

- **Resolution.** "Resolution-independent leak" if
  |Gamma(N192)-Gamma(N128)| <= 0.2 Gamma(N128) and
  |f_far(100), N192 - N128| <= 0.2 f_far(100, N128). Otherwise
  "resolution-dependent". N96 is reported alongside. If f_far(100)<1e-3 on
  both grids, the class is "no measurable leak".
- **m_TV dependence.** Gamma(m_TV=4)/Gamma(m_TV=1.5) <= 0.1 is "steep",
  >= 0.5 is "weak", otherwise "moderate".
- **Lifetime.** t_life/t_c for every run, with censored values given as lower
  bounds.
