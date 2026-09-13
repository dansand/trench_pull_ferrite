# Convergence & consistency — locked baseline (Tresca σ_Y=150 MPa, h=60 km, V=4 TN/m)

Trench deflection w_T, trench pull ΔGPE*, and the ΔN_D = ΔGPE* residual (identity consistency check).

## Spatial refinement (nsteps = 24)

| mesh nx × nz | w_T [m] | ΔGPE* [TN/m] | ΔN_D=ΔGPE* residual [%] |
|---|---|---|---|
| 400 × 24 | 3233 | 2.550 | 0.097 |
| 800 × 48 (baseline) | 3233 | 2.542 | 0.019 |
| 1200 × 72 | 3233 | 2.540 | 0.010 |

## Load-increment refinement (800 × 48)

| nsteps | w_T [m] | ΔGPE* [TN/m] | ΔN_D=ΔGPE* residual [%] |
|---|---|---|---|
| 12 | 3233 | 2.543 | 0.021 |
| 24 (baseline) | 3233 | 2.542 | 0.019 |
| 48 | 3233 | 2.542 | 0.019 |
