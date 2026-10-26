# Matched empty-pocket contraction control

Contract: [EMPTY_STIFF_CONTROL.md](../EMPTY_STIFF_CONTROL.md).
Raw compiled field run and independent radial history are in
`outputs/pockets/empty-stiff-control-01`.

The initial dilaton field and momentum match the N=256 gas run exactly.
Only the axion field/momentum are removed. Other parameters remain
L=64, dx=.25, dt=.0125, lambda=8, DeltaV=.01, m_TV=2, p=4,
R0=20, w=1, static a=1. The validated four-thread OpenMP build is used.

All predeclared early-control gates pass:

- Exact matched dilaton data and exact driver transfer.
- Maximum volume-radius error against the independent radial field PDE:
  .04435, below .25, through the trusted t=9.5 window (R_reference>=15).
- Maximum relative energy drift in that window: 9.86e-6, below 5e-4.
- All saved wall contour segments in the trusted window are measurable
  and at least four cells thick.

The radial reference first crosses R=15 at the sampled t=9.75. At that
time the matched gas pocket has R=19.469. This is controlled evidence that
gas changes the contraction for these stiff-wall parameters, not merely a
comparison to an older potential. The gas case remains above R=18.736 to
t=80, although its own candidate gates still fail as reported separately.

The empty run completes t=20 with recorded R=5.465. That later contraction
is outside the declared trusted-resolution window and is not used as an
accurate collapse-time result. None of these short runs establishes the
unchanged ten-crossing-time survival goal.
