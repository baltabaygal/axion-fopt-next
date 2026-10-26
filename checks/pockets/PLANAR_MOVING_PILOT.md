# First moving-wall coupled-field pilot

2026-10-18. A native CosmoLattice NDim=1 packet/matched-empty-wall pair
ran on one thread in 4.35 seconds, below a declared 15-second budget.
Evidence: `outputs/pockets/planar-moving-pilot-01/report.json`, raw inputs,
histories, final fields, matched wall output and frozen build identity.
No 3D experiment ran.

N=4096, L=512, dx=.125, dt=.0125, Minkowski time 0..180. Degenerate
lambda=8, mFV=.1, mTV=1, kappa=log(100), fa=10. The original positive-k
weak packet has amplitude .001, central omega=.5 and envelope width 16.
Both dilaton walls initially approach the false phase at speed .3 with
Lorentz-contracted profiles and canonical momentum -v*phi_x. Both fields
then evolve dynamically. The fixed diagnostic window is centre=88,
inner=72, outer=104; its geometry was checked against the reference first.

## Independent moving reference

`moving_packet_reference.py` boosts the independently solved static scattering
modes, rather than simply boosting an energy fraction. For target initial
position X0=256 and approaching speed u, use

xi=gamma(x-X0+u*t), tau=gamma(t+u*(x-X0)),
Omega=gamma(omega+u*k), K=gamma(k+u*omega).

For each mode a=Re[C psi(xi) exp(-i Omega tau)], the momentum is
Re[C gamma(u psi' - i Omega psi) exp(-i Omega tau)]. The spatial derivative
uses gamma(psi' - i u Omega psi). Exact free-wave boost tests check incident
lab phase, canonical momentum, gradient and reflected Doppler frequency.

The initial 61-mode sampling cannot be used unchanged as quadrature for a
single moving reflected packet: Doppler compression decreases its spatial
repeat distance and introduces an unwanted initial reflected copy. This was
caught before any compiled moving run. `planar-moving-reference-01` is
therefore an invalid initial-packet reference despite its original incomplete
checks. Fourfold quadrature (`reference-02`) removes the large alias but
fails the new, unchanged 1e-6 initial-core error gate. Eightfold quadrature
(`reference-03`) passes all gates. These earlier artifacts are retained;
no failed check is changed to passing. The runner now explicitly requires
initial-energy and initial-core checks, rejecting the earlier incomplete
reference even if its old aggregate flag is true.

The accepted reference's initial total energy differs from a single incident
packet by 1.82e-8 relatively; the initial core field/momentum differs from
the prepared input by at most 7.29e-7 in fixed incident scales. Initial
window coverage exceeds 99.9%; predicted final reflected-packet coverage is
99.92%. Reference evidence: `outputs/pockets/planar-moving-reference-03`.
This is a fixed-speed weak-field reference, not a predicted nonlinear wall
trajectory.

Lab energy must include the spectral Jacobian. The reflected wavenumber has
dk_out/dk=omega_out/omega. Thus modal returned energy is incident energy times
R_wall*(omega_out/omega), and the finite-band lab gain is 1.84392143. A
wall-frame R alone would miss this energy transfer. The central reflected
frequency is .921911, whereas the energy-weighted mean over the full packet
is .936026; these are distinct quantities.

## Actual coupled-field result

| Observable | Coupled fields | Reference |
|---|---:|---:|
| Final lab axion energy / initial | 1.84224584 | 1.84392143 |
| Energy-weighted reflected frequency | .93558725 | .93602631 |

Energy gain differs by .0909% and mean frequency by .0469%, passing the
predeclared 1% gates. The measured frequency uses the final Fourier energy
spectrum with the false-vacuum dispersion; it is dominated by the returned
packet. This does not measure a tiny transmitted component separately.

The matched packet-minus-empty change in dilaton energy is -.84225342 of
initial incident axion energy. Its cancellation with the axion gain leaves
7.58e-6, within the declared .002 matched-sector closure gate. The gain is
energy taken from a dynamical approaching wall, not energy creation.
The window encloses 99.913% of final axion energy. Its maximum integrated
axion continuity residual is 5.17e-5 of initial incident energy, passing
.002. The exact-zero control remains zero.

The final target wall is at x=202.095218, versus 202.093199 without axions.
Positive matched displacement is consistent with axion pressure opposing
contraction. Its contour-velocity estimate has interpolation bias and is
not treated as a precision material-velocity measurement.

All declared pilot gates pass, but the direct maximum waveform discrepancy
is .2251 in fixed initial field/momentum scales. This phase-sensitive error
is not hidden by the accurate integrated energy and frequency. The pilot
report retains `moving_scattering_validated: false`. Independent time,
space, amplitude and second-wall controls are required before promoting it
to a controlled moving benchmark. Neither this result nor the stationary
benchmark establishes long-lived 3D pockets, equilibrium scaling, expansion
or collision formation.
