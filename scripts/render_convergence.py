"""render_convergence.py — SI convergence + consistency table for the locked baseline (σ_Y=150, h=60, V=4).

Reads the `bench_*` models produced by `julia paper_models.jl convergence` and reports, per configuration:
  - w_T      trench deflection [m]
  - ΔGPE*    trench pull [TN/m]
  - residual the ΔN_D = ΔGPE* consistency residual [%]  (folds the identity check into the convergence table)

Two refinement axes: spatial (mesh) and load-increment count. Writes a Markdown table to
`tables/convergence.md` (+ `tables/convergence.csv`) and prints it.

    python scripts/render_convergence.py
"""
import sys; sys.path.insert(0, "analysis")
import os
import numpy as np
from gpe_analysis import Model, trench_pull, trench_deflection, write_table


def main():
    SB = "data/convergence"

    # The convergence data is regeneration-only (not shipped).  Refuse to run — and above all refuse to overwrite the
    # recorded table tables/convergence.md with an empty one — unless every bench_* model is present.
    _needed = ["bench_400x24", "bench_800x48", "bench_1200x72", "bench_800x48_ns12", "bench_800x48_ns48"]
    _missing = [n for n in _needed if not os.path.isfile(os.path.join(SB, n, "gpe_model.vtu"))]
    if _missing:
        raise SystemExit(f"render_convergence: {SB}/ is missing {len(_missing)} of {len(_needed)} bench_* models "
                         f"({', '.join(_missing)}).\nThe convergence data is not shipped; the recorded table is tables/convergence.md.\n"
                         "To regenerate:  julia --project=. model/paper_models.jl convergence   (then rerun this script).  Nothing was written.")

    rows_csv = []

    def row(name, label):
        d = os.path.join(SB, name)
        if not os.path.isfile(os.path.join(d, "gpe_model.vtu")):
            return f"| {label} | — | — | — (missing) |"
        m = Model(d)
        wT = trench_deflection(m)                                  # w_T on the trench column
        dG, dNd, _ = trench_pull(m)
        resid = abs(dG - dNd) / abs(dG) * 100
        rows_csv.append((name, label, wT, dG / 1e12, resid))
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

    out = "tables/convergence.md"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    write_table("convergence", ["model", "configuration", "w_T_m", "dGPE_TN", "identity_residual_pct"], rows_csv,
                script="scripts/render_convergence.py", models=[f"{SB}/{n}" for n, _ in spatial + loadinc])
    print("\n".join(lines))
    print("\nwrote", out, "and tables/convergence.csv")


if __name__ == "__main__":
    main()
