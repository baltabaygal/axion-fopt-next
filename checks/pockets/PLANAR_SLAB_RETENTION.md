# Wall-following slab retention: saved-endpoint audit

2026-10-19. The distributed slab's small fixed-window flux does not imply
negligible energy outside the pocket. An independent endpoint window tied
to both actual interfaces places about 3.1% of final axion energy outside
the slab window. No simulation was run for this audit.

The source fields and canonical units are those in PLANAR_DISTRIBUTED_PILOT.md.
The left and right phi=.5 interfaces and contour velocities are independently
located using original and reflected fields. Their mean gives window centre
and velocity; half their difference gives radius and radius velocity. A
quintic window transitions from one at R-padding to zero at R+padding.
The weight derivative includes both centre translation and radius change.
The fourth-order signed-bond Hamiltonian/current supplies energy, lab flux,
advection and instantaneous net outward flux. This is the same independent
continuity algebra previously tested against the native passive window.

| Endpoint window padding | Outside fraction of current axion energy |
|---:|---:|
| 2 | .0305582 |
| 4 | .0311825 |
| 8 | .0354306 |

At padding=4 the coarse-space and coarse-time values are .0311230 and
.0312034. Initial windows enclose the compact seed. These figures measure
weighted endpoint energy outside the interfaces, including penetration and
the transition region; they are not an asymptotic escaped fraction. The
signed local fourth-order gradient split also prevents interpreting every
individual cell as a positive particle energy. Padding sensitivity is
explicit, rather than selecting one boundary and calling it physical escape.

This resolves a concrete ambiguity: the original stationary window extends
into TV outside the shrinking slab, and can retain outgoing energy in that
buffer. Its net flux 2.64e-5 of initial energy is not inconsistent with the
~3% outside-interface endpoint energy. Reabsorption, transmission and a
lifetime cannot be distinguished from two snapshots. Full time-resolved
wall-following exchange/advection/flux must be measured before survival.

Source slab_retention.py; evidence
outputs/pockets/planar-distributed-endpoint-retention-01/report.json includes
three source report hashes, initial/final measures and padding sensitivity.
Two tests verify weight time derivatives against independent finite
differences for translation, contraction and combined motion, and verify
static-window limits. No heavy compute was used. Fine temporal/spatial and
outer-TV-buffer controls remain outstanding, along with the full 3D,
mass-radius scaling, FRW and formation objectives.
