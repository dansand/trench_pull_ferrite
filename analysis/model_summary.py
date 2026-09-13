"""model_summary.py — one row per production model: the trench pull and its companions, as one table.

The same numbers the notebook's §8 prints and the figures carry, in one place: plate thickness, applied load and
background N_D (from data/DATA_MANIFEST.json), trench deflection, the three columns, the pull ΔGPE*, ΔN_D, the
equilibrium identity residual, and the effective arm ΔGPE*/(Δρ g w_T).  Written to tables/model_summary.csv by
`reproduce.sh`; START_HERE §8 reads it.

    python analysis/model_summary.py
"""
import glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gpe_analysis import Model, reference_lines, trench_ref_km, trench_pull, write_table

SUITES = ["suite1_strength", "suite2_load", "suite3_background", "suite4_thickness"]
RHO_A, RHO_W, G = 3300.0, 1000.0, 9.81
DRHOG = (RHO_A - RHO_W) * G


def main():
    with open("data/DATA_MANIFEST.json") as f:
        manifest = json.load(f)["models"]
    rows = []
    for suite in SUITES:
        for d in sorted(glob.glob(f"data/{suite}/*/")):
            d = d.rstrip("/")
            if not os.path.isfile(os.path.join(d, "gpe_model.vtu")):
                continue
            m = Model(d); R = reference_lines(m); x_t = trench_ref_km(m)
            dG, dN, x_i = trench_pull(m); w_t = float(np.nanmax(m.topography()))
            par = manifest.get(d.replace("data/", ""), {}).get("parameters", {})
            eps_p = float(np.nanmax(m.array("plastic_strain"))) if m.has_field("plastic_strain") else 0.0
            H = float(par.get("hardening_H_Pa", 0.0))
            rows.append(dict(suite=suite, model=os.path.basename(d), strength=par.get("strength", ""), h_km=m.H / 1e3,
                             V_TN=par.get("V_TN_tuned", par.get("V_TN", float("nan"))), N_mem_TN=par.get("N_mem_TN", 0.0),
                             w_T_m=w_t, x_T_km=x_t, x_I_km=x_i, x_M_km=R["moment_max"], dGPE_TN=dG / 1e12, dND_TN=dN / 1e12,
                             identity_pct=abs(dG - dN) / abs(dG) * 100, arm_km=dG / (DRHOG * w_t) / 1e3, arm_over_h=dG / (DRHOG * w_t) / m.H,
                             max_plastic_strain=eps_p, hardening_H_Pa=H, hardening_increment_MPa=H * eps_p / 1e6))
    fields = ["suite", "model", "strength", "h_km", "V_TN", "N_mem_TN", "w_T_m", "x_T_km", "x_I_km", "x_M_km", "dGPE_TN", "dND_TN", "identity_pct", "arm_km", "arm_over_h", "max_plastic_strain", "hardening_H_Pa", "hardening_increment_MPa"]
    path = write_table("model_summary", fields, rows, script="analysis/model_summary.py", models=[f"data/{r['suite']}/{r['model']}" for r in rows])
    print(f"{'suite':18s} {'model':30s} {'h':>3s} {'V':>5s} {'N_mem':>5s} {'w_T[m]':>7s} {'ΔGPE*':>7s} {'ΔN_D':>7s} {'ident%':>7s} {'arm':>5s}")
    for r in rows:
        print(f"{r['suite']:18s} {r['model']:30s} {r['h_km']:3.0f} {r['V_TN']:5.2f} {r['N_mem_TN']:5.1f} {r['w_T_m']:7.0f} {r['dGPE_TN']:7.3f} {r['dND_TN']:7.3f} {r['identity_pct']:7.3f} {r['arm_km']:5.1f}")
    print("wrote", path)


if __name__ == "__main__":
    main()
