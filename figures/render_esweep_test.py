"""render_esweep_test.py — EXPLORATORY (TBD): does an elastic Young's-modulus sweep land on the
ΔGPE*-vs-w collapse line?  E is a FIXED-h (fixed-arm) knob for the deflection w: at fixed load V=4,
softer plates (low E) deflect more, stiffer less, but the moment arm z_np≈h/2 is unchanged.  So if the
trench pull is topographic — ΔGPE* ≈ Δρg·(h/2)·w — the E-sweep should slide ALONG the single line the
V-sweep and N_D-sweep already trace, NOT off it.  Writes data/reference_figures/esweep.png (Appendix B).

    /opt/anaconda3/envs/pyvista-env/bin/python render_esweep_test.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from gpe_analysis import Model, trench_pull
from render_gpe_correlation import SUITE2, SUITE3, SUITE1, gather, gather1, r2

# E-sweep at V=4, N_D=0: (E/E0 multiplier, path).  E0 = 70 GPa; baseline E×1 = elastic_deep_60km_V4 (Suite 1).
ES = "data/appendixB_esweep"
ESWEEP = [(0.125, f"{ES}/elastic_deep_60km_V4_E0p125"), (0.25, f"{ES}/elastic_deep_60km_V4_E0p25"),
          (0.5, f"{ES}/elastic_deep_60km_V4_E0p5"), (1.0, "data/suite1_strength/elastic_deep_60km_V4"),
          (2.0, f"{ES}/elastic_deep_60km_V4_E2"), (4.0, f"{ES}/elastic_deep_60km_V4_E4"),
          (8.0, f"{ES}/elastic_deep_60km_V4_E8")]
OUT = "data/reference_figures/esweep.png"
E0_GPA = 70.0


def gatherE(rows):
    mul, w, G = [], [], []
    for m_, p in rows:
        m = Model(p); mul.append(m_)
        w.append(m.topography().max() / 1e3); G.append(trench_pull(m)[0] / 1e12)
    return np.array(mul), np.array(w), np.array(G)


def main():
    V2, _, w2, G2 = gather(SUITE2)
    V3, N3, w3, G3 = gather(SUITE3)
    _, w1, G1 = gather1(SUITE1)
    mul, wE, GE = gatherE(ESWEEP)
    PURPLE = "#7b3fbf"

    # the EXISTING collapse line (V + N_D + rheology, fixed h) — the reference the E-sweep is tested against
    wref = np.concatenate([w2, w3, w1]); Gref = np.concatenate([G2, G3, G1])
    Rref, a, b = r2(wref, Gref)
    # how well does the E-sweep sit on THAT line?  (residual R² about the reference fit, not a refit)
    resid = GE - (a * wE + b)
    ss = 1 - np.sum(resid ** 2) / np.sum((GE - GE.mean()) ** 2)

    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.8), constrained_layout=True)

    # (A) ΔGPE* vs V — the E-sweep is a VERTICAL scatter at V=4 (same load, different w) spearing the V-line
    RV, av, bv = r2(V2, G2)
    ax[0].plot(V2, G2, "o-", color="#1f77b4", label="Suite 2 (vary V)")
    ax[0].plot([V2.min(), V2.max()], [av * V2.min() + bv, av * V2.max() + bv], "k--", lw=1, alpha=0.5)
    ax[0].scatter(V3, G3, c="#d62728", marker="s", zorder=5, label="Suite 3 (vary N$_D$)")
    ax[0].scatter(np.full_like(GE, 4.0), GE, c=PURPLE, marker="D", s=42, zorder=6, edgecolor="k",
                  linewidth=0.4, label="E-sweep (vary E, V=4)")
    ax[0].set_xlabel(r"applied load  V  [TN m$^{-1}$]"); ax[0].set_ylabel(r"$\Delta$GPE$^{*}$  [TN m$^{-1}$]")
    ax[0].set_title(f"(A)  vs load — E-sweep is a vertical spike at V=4", fontsize=10); ax[0].grid(alpha=0.25)
    ax[0].legend(fontsize=8.5)

    # (B) ΔGPE* vs w — does the E-sweep land ON the collapse line?
    ax[1].plot(w2, G2, "o", color="#1f77b4", alpha=0.8, label="Suite 2 (V-driven w)")
    ax[1].scatter(w3, G3, c="#d62728", marker="s", zorder=4, alpha=0.85, label="Suite 3 (N$_D$-driven w)")
    ax[1].scatter(w1, G1, c="#2ca02c", marker="^", s=48, zorder=4, edgecolor="k", linewidth=0.4,
                  label="Suite 1 (rheology)")
    xw = np.array([min(w2.min(), wE.min()), max(w2.max(), wE.max())])
    ax[1].plot(xw, a * xw + b, "k--", lw=1, alpha=0.7, label=f"collapse fit (V+N$_D$+rheo)  R$^2$={Rref:.3f}")
    ax[1].plot(wE, GE, "-", color=PURPLE, lw=1, alpha=0.5, zorder=5)
    ax[1].scatter(wE, GE, c=PURPLE, marker="D", s=46, zorder=6, edgecolor="k", linewidth=0.4,
                  label=f"E-sweep  (on-line R$^2$={ss:.3f})")
    for m_, x, y in zip(mul, wE, GE):        # annotate E multiplier
        ax[1].annotate(f"{m_:g}×", (x, y), fontsize=6.5, color=PURPLE, xytext=(4, 4), textcoords="offset points")
    ax[1].set_xlabel(r"trench deflection  w  [km]"); ax[1].set_ylabel(r"$\Delta$GPE$^{*}$  [TN m$^{-1}$]")
    ax[1].set_title(f"(B)  vs topography — does the E-sweep collapse too?", fontsize=10); ax[1].grid(alpha=0.25)
    ax[1].legend(fontsize=8)

    fig.suptitle(r"E-sweep test (V=4, N$_D$=0):  varying Young's modulus slides $w$ at fixed arm — "
                 "does $\\Delta$GPE$^{*}$ stay on the collapse line?", fontsize=11)
    fig.savefig(OUT, dpi=140); print("wrote", OUT)
    slope_geo = 2300 * 9.81 * 30000 / 1e12 * 1e3          # Δρg·(h/2) in TN/m per km
    print(f"  collapse-line slope (fit):     {a:.3f} TN/m per km   (analytic Δρg·h/2 = {slope_geo:.3f})")
    print(f"  E-sweep on-line R^2:           {ss:.4f}")
    print("  E×mult   E[GPa]    w[km]   ΔGPE*[TN/m]   line-pred   resid")
    for m_, x, y, r in zip(mul, wE, GE, resid):
        print(f"  {m_:6g}  {m_*E0_GPA:7.1f}  {x:6.2f}   {y:9.3f}   {a*x+b:8.3f}  {r:+.3f}")


if __name__ == "__main__":
    main()
