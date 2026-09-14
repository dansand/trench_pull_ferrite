"""edge_exclusion.py — the loaded-edge sensitivity the SI quotes (Text S1, model formulation): what excluding the
first ~10 km changes in the deflection-normalised effective arm.

Restored (as a table, no figure) from the retired render_edge_si.py.  Per Suite-1 model: the arm ΔGPE*(x)/(Δρ g w(x))
profile from the trench inboard, on deformed columns; the edge window L_e where the top-10-km equivalent density has
decayed to its inboard level; the plateau arm (median over [L_e, 2L_e]); and the change from the raw trench arm to the
plateau.  Writes tables/edge_exclusion.csv.

    python analysis/edge_exclusion.py
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gpe_analysis import Model, trench_ref_km, deformed_shear_gradient, reference_lines, write_table, RHO_M, RHO_W, G

DR, THR = RHO_M - RHO_W, 0.15
S1 = "data/suite1_strength/"
SUITE = [("elastic", S1 + "elastic_deep_60km_V4"), ("Tresca", S1 + "tresca_deep_150_60km_V4"),
         ("DD-VM asym", S1 + "dd_vm_asym_60km_V4"), ("DD-VM sym", S1 + "dd_vm_sym_60km_V4")]


def profile(p):
    m = Model(p); xt = trench_ref_km(m)
    xi = reference_lines(m)["isostatic"]; SzzI = m.deformed_resultants(xi)[0]; top = m.topography()
    XS = np.linspace(max(xt, 0.4), 90, 200); arm, skin = [], []
    for x in XS:
        L, tau = deformed_shear_gradient(m, x); d = (L["z"] - L["z"][0]) / 1e3; skin.append(np.mean(tau[d <= 8.0]) / G)
        Szz = m.deformed_resultants(x)[0]; w = np.interp(x, m.xkm, top)
        arm.append((Szz - SzzI) / (G * DR * w) / 1e3 if w > 50 else np.nan)         # direct arm ΔGPE*(x)/(Δρ g w) [km]
    arm, skin = np.array(arm), np.array(skin)
    sfar = np.median(skin[(XS >= 40) & (XS <= 80)])
    dec = np.where((XS > xt) & (np.abs(skin - sfar) < THR * abs(skin[0] - sfar)))[0]
    Le = float(XS[dec[0]]) if dec.size else 10.0
    win = (XS >= Le) & (XS <= 2 * Le)
    plat = float(np.nanmedian(arm[win]))
    return dict(h_km=m.H / 1e3, raw_trench_arm_km=float(arm[0]), edge_window_km=Le, plateau_arm_km=plat,
                plateau_min_km=float(np.nanmin(arm[win])), plateau_max_km=float(np.nanmax(arm[win])),
                change_raw_to_plateau_pct=100 * (plat - arm[0]) / arm[0], plateau_arm_over_h=plat / (m.H / 1e3),
                arm_at_40km_km=float(np.interp(40, XS, arm)), arm_at_60km_km=float(np.interp(60, XS, arm)))


def main():
    rows = []
    for lab, p in SUITE:
        r = profile(p); r["label"] = lab; r["model"] = os.path.basename(p); rows.append(r)
        print(f"{lab:10s} raw arm {r['raw_trench_arm_km']:.1f} km  edge window {r['edge_window_km']:.1f} km  plateau {r['plateau_arm_km']:.1f} km "
              f"({r['plateau_arm_over_h']:.2f} h)  change {r['change_raw_to_plateau_pct']:+.1f} %")
    fields = ["label", "model", "h_km", "raw_trench_arm_km", "edge_window_km", "plateau_arm_km", "plateau_min_km", "plateau_max_km",
              "change_raw_to_plateau_pct", "plateau_arm_over_h", "arm_at_40km_km", "arm_at_60km_km"]
    print("wrote", write_table("edge_exclusion", fields, rows, script="analysis/edge_exclusion.py", models=[p for _, p in SUITE]))


if __name__ == "__main__":
    main()
