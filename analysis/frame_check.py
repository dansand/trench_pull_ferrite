"""frame_check.py — the finite-rotation frame check the SI quotes (Text S2, 'A final note concerns the frame').

N_D is computed on vertical columns and is exact at any slope; identifying it with a plate-aligned resultant is a
small-slope approximation, N_D^plate = N_D cos 2θ + 2V sin 2θ with θ the local surface slope and V = ∫τzx dz.  Per
production model this script measures: the maximum surface slope; the maximum departure of cos 2θ from 1; the mixing
term 2V sin 2θ at the trench column; the distance inboard at which it first falls below 0.02 TN/m and beyond which it stays there; and the slope and mixing
term at the first isostatic column.  V is the reference-grid Cauchy shear resultant (Model.V, the column-locating
estimator), θ from the deformed surface.  Writes tables/frame_check.csv.

    python analysis/frame_check.py
"""
import glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gpe_analysis import Model, reference_lines, trench_ref_km, write_table

SUITES = ["suite1_strength", "suite2_load", "suite3_background", "suite4_thickness"]
THRESH_TN = 0.02


def measure(mdir):
    m = Model(mdir)
    xs = m.x + m.array("u", 0)[:, 0]; w = m.topography()                     # deformed surface x [m], deflection [m]
    theta = np.arctan(np.gradient(w, xs))                                     # local surface slope [rad]
    V = m.V()                                                                 # ∫τzx dz [N/m] on the reference grid
    mix = 2 * V * np.sin(2 * theta) / 1e12                                    # 2V sin2θ [TN/m]
    xk = m.xkm; x_t = trench_ref_km(m); x_i = reference_lines(m)["isostatic"]
    it = np.argmin(np.abs(xk - x_t)); ii = np.argmin(np.abs(xk - x_i))
    beyond = np.where((xk > x_t) & (np.abs(mix) < THRESH_TN))[0]
    first = xk[beyond[0]] if beyond.size else float("nan")                 # first x inboard where |2V sin2θ| < 0.02
    # first x from which the term STAYS below the threshold
    stay = next((xk[j] for j in beyond if np.all(np.abs(mix[j:]) < THRESH_TN)), float("nan"))
    return dict(model=os.path.basename(mdir), h_km=m.H / 1e3, max_slope_deg=float(np.degrees(np.nanmax(np.abs(theta)))),
                max_cos2theta_departure_pct=float(100 * np.nanmax(1 - np.cos(2 * theta))),
                mixing_term_at_trench_TN=float(mix[it]), mixing_first_below_0p02_km=float(first - x_t), mixing_below_0p02_beyond_km=float(stay - x_t),
                slope_at_xI_deg=float(np.degrees(abs(theta[ii]))), mixing_term_at_xI_TN=float(mix[ii]), x_I_km=x_i)


def main():
    rows = []
    for suite in SUITES:
        for d in sorted(glob.glob(f"data/{suite}/*/")):
            if os.path.isfile(os.path.join(d, "gpe_model.vtu")):
                r = measure(d.rstrip("/")); r["suite"] = suite; rows.append(r)
    fields = ["suite", "model", "h_km", "max_slope_deg", "max_cos2theta_departure_pct", "mixing_term_at_trench_TN",
              "mixing_first_below_0p02_km", "mixing_below_0p02_beyond_km", "slope_at_xI_deg", "mixing_term_at_xI_TN", "x_I_km"]
    path = write_table("frame_check", fields, rows, script="analysis/frame_check.py", models=[f"data/{r['suite']}/{r['model']}" for r in rows])
    for r in rows:
        print(f"{r['model']:32s} slope max {r['max_slope_deg']:.2f}°  cos2θ dep {r['max_cos2theta_departure_pct']:.2f}%  mix@trench {r['mixing_term_at_trench_TN']:+.3f}  "
              f"<0.02 beyond {r['mixing_below_0p02_beyond_km']:.0f} km  slope@xI {r['slope_at_xI_deg']:.2f}°  mix@xI {r['mixing_term_at_xI_TN']:+.3f}")
    print("wrote", path)


if __name__ == "__main__":
    main()
