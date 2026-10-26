# Matched empty-pocket control for the stiff-wall state

Declared before run, 2026-10-17. Match the fine stiff-wall case's initial
dilaton field and momentum exactly, remove only the axion field/momentum,
and compare the early contraction to an independent radial field PDE.
This checks the role of gas for these parameters rather than borrowing
the older lambda=2 empty control.

N=256,L=64,dx=.25,dt=.0125,t_end=20,lambda=8,DeltaV=.01,m_TV=2,
p=4,R0=20,w=1,a=1. Histories every.25, slices every1, full states every10.
Use the validated accurate OpenMP executable, four threads. No expansion
or absorber. The axion is identically zero, so the phi^4 coupling is inert.

Reference radial PDE: r_max=32,dr=.025,dt=.002, samples every.25,
same tanh profile, zero initial momentum, phi=1 at the outer boundary.

Require exact initial dilaton/momentum agreement with the gas run, exact
driver transfer, energy drift<5e-4 and radius agreement<.25 while the
reference radius remains >=15. This early window is chosen because a
relativistically collapsing empty wall contracts far below the initial
width; do not use its unresolved implosion as an accurate collapse time.
Report widths at every saved plane; require >=4 cells and unambiguous
walls in the trusted window. Outside it, record failures without claiming
validated late collapse. Report the earliest sampled crossing of R=15
and compare the gas radius at the same time.
