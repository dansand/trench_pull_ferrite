"""render_mp_benchmark.py — SI plasticity benchmark: elasto-plastic moment–curvature to the fully-plastic M_p.

The idealized clamped–guided beam (`idealized_beam.jl`: displacement-controlled, no foundation, σ_zz ≈ 0,
Tresca, perfect plasticity H=0) has a moment that varies along its length, so ONE final state sweeps every
section from elastic (M≈0, midspan) through first yield to nearly fully plastic (near the ends). The
elastic-perfectly-plastic M(κ) curve is UNIVERSAL in normalized coordinates, so normalizing each section by
its LOCAL yield stress (raised by the end-protection ramp) collapses clean and ramped sections onto one curve:

    M/M_p  vs  κ/κ_y ,   M_p = σ_Y h²/4 ,  M_y = σ_Y h²/6 = (2/3)M_p ,  κ_y = 2σ_Y/(E'h)
    analytic:  M/M_p = (2/3)(κ/κ_y)                for κ < κ_y   (elastic)
               M/M_p = 1 − (1/3)(κ_y/κ)²           for κ ≥ κ_y   (elastic-perfectly-plastic)

  (a) M(κ): FE sections vs the universal analytic curve   (b) σ_xx(z) truncating at ±σ_Y (the M_p mechanism)

STRESS FRAME — material (2nd Piola–Kirchhoff S), BY DESIGN.  M and κ are section quantities of the material beam
(reference thickness h, fibre stress on the reference section), so the exported `sigma_xx` is read directly.  The
CLAUDE.md "always Cauchy" rule is for the deformed-column production pipeline, not for this large-rotation
material-section benchmark; pushing S to Cauchy mixes axial and shear under the rotation (2026-09-11 regression:
mean/max misfit 0.09/0.85 % → 0.57/5.24 %).  Prints the SI Text S1 numbers and ASSERTS them.

    python scripts/render_mp_benchmark.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from matplotlib import cm
from gpe_analysis import Model, write_table


def main():
    E, NU = 60.0, 0.25; EP = E / (1 - NU**2)                # match idealized_beam.jl
    m = Model("data/idealized_beam")
    z = m.z; zmid = z.mean(); H = z.max() - z.min(); Iz = H**3 / 12
    x = m.x; L = x.max()
    sxx = m.array("sigma_xx [Pa]"); yld = m.array("yielded")     # material-frame S, see the docstring (NOT cauchy_fields)
    sY = 1.0 * (1.0 + 8.0 * (np.exp(-x / (0.05 * L)) + np.exp(-(L - x) / (0.05 * L))))   # the end-protection ramp
    yf = np.nanmean(yld > 0.5, axis=1)                      # through-thickness yield fraction per section

    M = np.array([np.trapz(sxx[i] * (z - zmid), z) for i in range(m.Nx)])                # bending moment
    kap = np.full(m.Nx, np.nan)
    for i in range(m.Nx):                                   # curvature from the elastic-core slope of σxx(z)
        s = sxx[i]; core = np.abs(s) < 0.9 * sY[i]
        if core.sum() >= 4:
            kap[i] = np.polyfit(z[core] - zmid, s[core], 1)[0] / EP

    M_p = sY * H**2 / 4; kap_y = 2 * sY / (EP * H)
    Mn = np.abs(M) / M_p; kn = np.abs(kap) / kap_y
    analytic = lambda k: np.where(k < 1, (2 / 3) * k, 1 - (1 / 3) / np.maximum(k, 1e-9)**2)

    ok = np.isfinite(kn) & (kn > 0.15) & (kn < 8) & (np.abs(M) > 0.02)
    resid = np.abs(Mn[ok] - analytic(kn[ok]))
    kmax = kn[ok].max()
    print(f"sections on the curve: {ok.sum()}   mean |M/Mp − analytic| = {resid.mean()*100:.2f}%   max {resid.max()*100:.2f}%")
    print(f"max κ/κ_y = {kmax:.2f}  →  M/M_p = {Mn[ok][np.argmax(kn[ok])]:.3f}  (analytic {analytic(kmax):.3f})")

    # --- acceptance: the SI numbers, asserted (not just printed — the 2026-09-11 regression slipped past a print) ---
    assert resid.mean() <= 0.0015, f"mean M–κ misfit {resid.mean()*100:.2f}% > 0.15% — wrong stress frame?"
    assert resid.max() <= 0.01,    f"max M–κ misfit {resid.max()*100:.2f}% > 1%"
    assert ok.sum() == 148,        f"{ok.sum()} sections on the curve, SI states 148"
    print("benchmark assertions passed")
    write_table("benchmark_mp", ["quantity", "value", "unit"], [
        ("sections_on_curve", int(ok.sum()), "-"), ("mean_misfit_M_over_Mp", resid.mean() * 100, "%"),
        ("max_misfit_M_over_Mp", resid.max() * 100, "%"), ("kappa_max_over_kappa_y", kmax, "-"),
        ("M_over_Mp_at_kappa_max", float(Mn[ok][np.argmax(kn[ok])]), "-"), ("analytic_M_over_Mp_at_kappa_max", float(analytic(kmax)), "-")],
        script="scripts/render_mp_benchmark.py", figure="figures/benchmark_mp.png", models=["data/idealized_beam"])

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)

    # (a) moment–curvature
    ka = np.linspace(0, max(3.0, kmax * 1.05), 400)
    ax[0].plot(ka, analytic(ka), color="#b2182b", lw=3.0, alpha=0.5, label="analytic (elastic–perfectly-plastic)")
    sc = ax[0].scatter(kn[ok], Mn[ok], c=yf[ok], cmap="plasma", s=16, ec="k", lw=0.2, zorder=3, label="FE sections")
    ax[0].axhline(1.0, color="0.5", lw=0.9, ls=":"); ax[0].text(0.15, 1.01, "$M_p=\\sigma_Y h^2/4$", fontsize=8, color="0.4")
    ax[0].plot(1, 2/3, "o", ms=7, mfc="none", mec="k", mew=1.2)
    ax[0].annotate("first yield\n$M_y=\\frac{2}{3}M_p$", (1, 2/3), xytext=(1.15, 0.42), fontsize=8,
                   arrowprops=dict(arrowstyle="->", color="0.4"))
    ax[0].set_xlim(0, max(3.0, kmax * 1.05)); ax[0].set_ylim(0, 1.08)
    ax[0].set_xlabel(r"normalized curvature  $\kappa/\kappa_y$"); ax[0].set_ylabel(r"normalized moment  $M/M_p$")
    ax[0].set_title("(a)  moment–curvature to the plastic moment", fontsize=10.5)
    ax[0].annotate(f"mean dev {resid.mean()*100:.1f}%\nreaches {Mn[ok][np.argmax(kn[ok])]:.2f}$M_p$ at $\\kappa/\\kappa_y$={kmax:.1f}",
                   (0.05, 0.9), fontsize=8, color="0.3", va="top")
    ax[0].legend(fontsize=8.5, loc="lower right"); ax[0].grid(alpha=0.2)
    fig.colorbar(sc, ax=ax[0], shrink=0.8, label="yield fraction")

    # (b) σ_xx(z) truncating at ±σ_Y as yielding spreads
    zeta = (z - zmid) / (H / 2)
    targets = [0.4, 0.8, 1.2, 1.6, 2.0]
    cmap = cm.plasma; nrm = plt.Normalize(0, max(yf[ok].max(), 1e-9))
    for kt in targets:
        j = ok & (np.abs(kn - kt) < 0.12)
        if j.sum():
            i = np.where(j)[0][np.argmin(np.abs(kn[j] - kt))]
            ax[1].plot(sxx[i] / sY[i], zeta, color=cmap(nrm(yf[i])), lw=1.8)
    ax[1].axvline(1, color="0.6", lw=0.8, ls="--"); ax[1].axvline(-1, color="0.6", lw=0.8, ls="--")
    ax[1].text(1.03, -0.9, r"$\pm\sigma_Y$", fontsize=8, color="0.4"); ax[1].set_xlim(-1.2, 1.2); ax[1].set_ylim(1, -1)
    ax[1].set_xlabel(r"$\sigma_{xx}/\sigma_Y$"); ax[1].set_ylabel(r"normalized depth  $(z-h/2)/(h/2)$")
    ax[1].set_title("(b)  stress truncates at $\\pm\\sigma_Y$ (the $M_p$ mechanism)", fontsize=10.5)
    ax[1].grid(alpha=0.2)

    fig.suptitle("Plasticity benchmark: elasto-plastic bending vs the analytic moment–curvature curve", fontsize=11.5)
    OUT = "figures/benchmark_mp.png"
    fig.savefig(OUT, dpi=140); print("wrote", OUT)


if __name__ == "__main__":
    main()
