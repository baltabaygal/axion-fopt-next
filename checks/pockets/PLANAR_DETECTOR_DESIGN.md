# Detector redesign from the finite-band reference

2026-10-18. No lattice simulation ran. The original pilot's 98% detector
coverage gate remains failed and its outputs are unchanged.

`design_planar_detector.py` declares the new geometry and gates before
computing the continuum packet. The reference uses the same 61 positive-k
incident coefficients and independently integrated complex scattering modes.
It computes local continuum energy/current over t=0..220, with samples every
.5, in batches of eight times to keep memory small. It models a fixed
degenerate wall in the weak-field limit; it is not a coupled-field run.

Evidence: `outputs/pockets/planar-detector-design-01/report.json`, the declared
contract, original/proposed reference energy/flux curves, complex mode
coefficients and source hashes. The reference calculation took under one
second. Detector radii are measurement parameters only; no potential,
initial field, simulation grid, time window or physical success requirement
was changed.

| Quantity, normalized to reference incident energy | Original window | Proposed window |
|---|---:|---:|
| Centre / inner / outer | 96 / 80 / 112 | 104 / 88 / 120 |
| Initial packet energy enclosed | 0.99999980 | 0.9999999999 |
| Final packet energy enclosed | 0.93720067 | 0.99745006 |
| Minimum enclosed during encounter | 1.75e-7 | 8.72e-5 |
| Outward crossing | 0.99999963 | 0.99991270 |
| Returned crossing | 0.93717623 | 0.99736022 |

Predeclared new-reference requirements are initial coverage >99.9%, final
coverage >98%, minimum enclosed during encounter <1%, and relative energy
balance residual <2e-3. The proposed window passes all. Its largest integrated
energy residual is 3.011e-5. It allows essentially all the incident energy to
leave before return, so improved endpoint coverage is not achieved merely by
keeping the interaction permanently inside the detector.

The quoted fractions use the continuum reference's own incident energy,
not a discretely sampled pilot normalization. The earlier reference used
endpoint-shared discrete gradient energy, so tiny numerical differences
between those reference fractions are expected and explicitly distinguished.

## Saved-endpoint check, without retrospective success

Only after the reference design was selected, the existing compiled pilot's
final state was remeasured offline with the proposed window and the independent
semi-discrete diagnostic. It encloses 0.99717243 of that pilot's initial
axion energy. Evidence:
`outputs/pockets/planar-detector-endpoint-audit-01/report.json`, with hashes of
the old final state, old report, reference design and diagnostic implementation.

This is an endpoint-only check, not a new run. There is no new-window flux
history in the old artifacts, so no new reflection coefficient or pilot
success is claimed. The old detector coverage failure remains false in its
original report. The new geometry must be used prospectively.

`run_planar_packet_pilot.py` now optionally accepts `--detector-design` with
the passing reference report. It checks that the design's grid, box and
end time match the pilot and saves the design's file identity and hash in
its new contract. The default still retains the original geometry. The
runtime limit, 98% endpoint coverage, incoming-crossing and energy-accounting
gates are unchanged. The runner has not been launched with the new design.

Next run one bounded new-detector case, then separately refine time, space,
amplitude and boundaries. Even a passing detector pilot will not establish
a converged transmission coefficient or three-dimensional pocket survival.
