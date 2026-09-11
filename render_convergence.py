"""render_convergence.py — SI convergence + consistency table for the locked baseline (σ_Y=150, h=60, V=4).

Reads the `bench_*` models produced by `julia paper_models.jl convergence` and reports, per configuration:
  - w_T      trench deflection [m]
  - ΔGPE*    trench pull [TN/m]
  - residual the ΔN_D = ΔGPE* consistency residual [%]  (folds the identity check into the convergence table)

Two refinement axes: spatial (mesh) and load-increment count. Writes a Markdown table to
`out/paper/CONVERGENCE.md` and prints it.

    /opt/anaconda3/envs/pyvista-env/bin/python render_convergence.py
"""
import sys; sys.path.insert(0, ".")
import os
import numpy as np
from gpe_analysis import Model, trench_pull

SB = "out/paper/convergence"

def row(name, label):
    d = os.path.join(SB, name)
    if not os.path.isfile(os.path.join(d, "gpe_model.vtu")):
        return f"| {label} | — | — | — (missing) |"
    m = Model(d)
    wT = float(np.nanmax(m.topography()))
    dG, dNd, _ = trench_pull(m)
    resid = abs(dG - dNd) / abs(dG) * 100
    return f"| {label} | {wT:.0f} | {dG/1e12:.3f} | {resid:.3f} |"

spatial = [("bench_400x24",  "400 × 24"),
           ("bench_800x48",  "800 × 48 (baseline)"),
           ("bench_1200x72", "1200 × 72")]
loadinc = [("bench_800x48_ns12", "12"),
           ("bench_800x48",      "24 (baseline)"),
           ("bench_800x48_ns48", "48")]

lines = ["# Convergence & consistency — locked baseline (Tresca σ_Y=150 MPa, h=60 km, V=4 TN/m)", "",
         "Trench deflection w_T, trench pull ΔGPE*, and the ΔN_D = ΔGPE* residual (identity consistency check).",
         "", "## Spatial refinement (nsteps = 24)", "",
         "| mesh nx × nz | w_T [m] | ΔGPE* [TN/m] | ΔN_D=ΔGPE* residual [%] |",
         "|---|---|---|---|"]
lines += [row(n, l) for n, l in spatial]
lines += ["", "## Load-increment refinement (800 × 48)", "",
          "| nsteps | w_T [m] | ΔGPE* [TN/m] | ΔN_D=ΔGPE* residual [%] |",
          "|---|---|---|---|"]
lines += [row(n, l) for n, l in loadinc]

out = "out/paper/CONVERGENCE.md"
with open(out, "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print("\nwrote", out)
