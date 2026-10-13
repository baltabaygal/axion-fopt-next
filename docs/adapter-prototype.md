# Arrival-map C++ feasibility prototype

This is a compiled host-side transfer and potential evaluator, not a working
CosmoLattice model or a field-evolution benchmark. No upstream code is vendored.

## Run

```sh
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python examples/export_one_bubble.py
c++ -std=c++20 -Wall -Wextra -Werror adapters/cosmolattice/probe.cpp -o /tmp/afopt-probe
/tmp/afopt-probe outputs/one_bubble.afopt 2 .5 2 .4 .2 1.5
```

The probe prints program mass squared, potential derivative, and potential
energy for each cell. Arguments are map, program stage time, conformal origin,
omegaStar, theta, false mass, true mass. Tests compare C++ with the independent
Python API, including periodic one-bubble geometry, nonmonotonic time queries,
exact arrival, unconverted cells, and malformed files.

## Binary contract

`AFOPT001` (8 ASCII bytes), n (little-endian uint32), box length (little-endian
IEEE754 float64), then n^3 arrival times (same float64). Header length is 20 bytes.
Grid axes are x,y,z; C order means z varies fastest. Points lie at
(x+1/2,y+1/2,z+1/2)*L/n. Positive infinity means never converted. NaN and negative
infinity are invalid. The reader rejects missing or trailing bytes and bounds
allocation (default maximum 256^3 cells; explicitly adjustable in the C++ API).

The format describes geometry only; masses, cosmology, and solver units remain
run parameters. A JSON sidecar in the example records the deterministic catalogue
and checksum. The C++ probe does not validate this sidecar checksum.

## Unit restriction

For the initial adapter choose CosmoLattice alpha=1, fStar=f_axion, hence the
program scalar is theta and program elapsed time u=omegaStar*(tau-tau_origin).
The transferred box becomes L_program=omegaStar*L. At each force stage:

    tau = tau_origin + u_stage/omegaStar
    m_program^2 = (m_phase/omegaStar)^2
    V_program = m_program^2*(1-cos(theta))
    dV_program/dtheta = m_program^2*sin(theta)

The C++ evaluator applies no scale-factor power: CosmoLattice's scalar momentum
kernel supplies a^(3+alpha), and its field drift supplies a^(alpha-3). The
resulting second-order conformal equation has a^2 multiplying the force. This
mapping does not apply to arbitrary alpha without integrating the time map.

## Findings in pinned upstream source

Inspected official tag v2.0.0, commit
`0dad16ec2c894a7c6abbea0e98734126b3470b36`, on 2026-10-13.

- `models/lphi4.h` defines potentialTerms, potDeriv, potDeriv2 without time inputs.
- `include/CosmoInterface/evolvers/kernels/scalarsingletkernels.h` receives an
  EoM kernel type containing tMinust0, but calls Potential::derivS(model,n)
  without forwarding it.
- `include/CosmoInterface/evolvers/rk2nstorage.h` caches tMinust0 at the beginning
  and after each stage, without a stage-time increment. Stage field states and
  fixed-background scale factors advance separately. Reusing that cached time
  for a new nonautonomous potential would evaluate the wrong background.
- `include/CosmoInterface/evolvers/leapfrog.h` scalar kicks use the default kernel
  overload, whose EoM time is zero. Measurement synchronization also kicks.

Therefore a model-only mass update once per outer step is insufficient. Add an
explicit background-time hook at every force and measurement evaluation, using
stage abscissae derived from the selected integrator's actual stored-state
convention. Do not infer the abscissae from the coefficient name alone. Validate
that hook with delayed homogeneous turn-on before any stochastic run.

## Remaining gate

Compile a real TempLat-backed model; establish node/cell offset compatibility;
transfer into rank-local lattice storage; add correctly timed force, energy and
mass evaluators; test homogeneous ODE limits and wall-crossing convergence.
GPU/device expressions and MPI decomposition are not covered by this host-side
reader. Start with fixed expansion: external energy input cannot be treated as
closed-system self-consistent Friedmann evolution without a source budget.

Upstream source: https://github.com/cosmolattice/cosmolattice/tree/v2.0.0
