# Homogeneous PQ comparison

Accepted run: `outputs/pq_strings/homogeneous-03`. All 87 numerical gates passed.

The full complex field approaches the fixed-radius theta-only solution as the radial mass increases. Both runs use the same radiation expansion, constant physical axion mass, and finite-start initial angle and velocity. The Bessel function sets initial data; it is not used as the nonlinear evolution reference.

| Input amplitude | m_r / m_a | Max phase difference / amplitude | Max absolute radial excursion |
|---|---|---|---|
| 0.01 | 2 | 2.13196e-05 | 1.46339e-05 |
| 0.01 | 8 | 6.08157e-07 | 1.92793e-06 |
| 0.01 | 32 | 3.53614e-08 | 1.51284e-07 |
| 0.8 | 2 | 0.0757031 | 0.0856376 |
| 0.8 | 8 | 0.00370802 | 0.011921 |
| 0.8 | 32 | 0.000215411 | 0.000924099 |

Parameters: f_a=m_a=1, chi=1, a=tau=1+u, u=0..8, L=2*pi, N=4. Radial lambda=2,32,512 is constant. Initial r=1 and r_prime=0. Angle and conformal angular velocity reuse checks/homogeneous.py at the finite start. Momenta are pi_X,Y=a^2*(X_prime,Y_prime). Errors are normalized to the fixed input amplitude, including momentum errors; they are not relative errors at zero crossings.

At dt=0.00005 and 0.000025, every compiled trajectory agrees with an independent polar-coordinate ODE within the declared field and final-momentum tolerances. Field error ratios are approximately four. Initial/final local spreads vanish, the gradient energy vanishes, and an N=8 repeat passes grid invariance. Grid invariance of a homogeneous solution does not establish spatial convergence for strings.

Earlier raw trials homogeneous-01 and homogeneous-02 are retained as failed numerical checks: their coarser timesteps failed the unchanged momentum gate for the heaviest larger-angle case. No acceptance tolerances were relaxed. The final run reduced both timesteps.

The finite-angle agreement is approximate. Radial motion remains physical, and the observed discrepancy is not a solver failure. The small-angle case suppresses radial feedback much more strongly.

Reproduce with the command in ../../README.md. Raw inputs, binary states, solver logs, histories, comparison CSVs and provenance remain in the ignored output directory. Compact contract and acceptance details accompany this file.
