"""render_benchmark.py — SI methodology benchmark: the finite-strain PRODUCTION code vs the analytic
beam-on-elastic-foundation (Hetényi) solution, on one page.

Uses the committed elastic baseline `suite1_strength/elastic_deep_60km_V4` (E=70 GPa, ν=0.25, h=60 km, Winkler
k=Δρg, end shear P) — the same finite-strain formulation used for every production model, run with σ_Y→∞.

  (a) deflection   w(x')  FE vs  w = (2Pβ/k) e^{-βx'} cos βx'
  (b) shear force  Q(x')  FE vs  Q = −P e^{-βx'} (cos βx' − sin βx')
  (c) shear stress τ_zx(z) FE vs the parabola  1.5 (Q/h)[1 − ((z−h/2)/(h/2))²]

β = (k/4D)^{1/4}, D = E'h³/12, E' = E/(1−ν²).  At h/α ≈ 0.48 (α = 1/β) the FE deflection sits a ~1.7% ABOVE
Euler–Bernoulli — the expected thick-plate (transverse-shear) offset, NOT residual error — annotated as such.
Prints the numbers for SI Text S1 and ASSERTS them (shear parabola, V misfit, end resultant).

STRESS FRAME — material (2nd Piola–Kirchhoff S), BY DESIGN.  The Hetényi comparison and the shear parabola are
formulated on the material section: reference thickness h, resultants per unit reference length.  The exported
`sigma_*` fields are exactly that stress, so they are read directly.  The CLAUDE.md rule "always go through
Model.cauchy_fields()" governs the deformed-column PRODUCTION pipeline (trench pull, ΔGPE*, ρ̂), where the lab-axis
Cauchy stress on the deformed geometry is what the force balance needs; it does NOT apply to these material-section
benchmarks.  Pushing S forward to Cauchy here mixes axial and shear components under the rotation and breaks the
comparison (seen 2026-09-11: parabola misfit 0.39 % → 346 %, end resultant 4.000 → 3.950 TN/m).

    python scripts/render_benchmark.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from gpe_analysis import Model, trench_ref_km, write_table


def main():
    # --- config constants (match paper_models.jl) ---
    E, NU, H = 70.0e9, 0.25, 60.0e3
    RHO_A, RHO_W, G = 3300.0, 1000.0, 9.81
    K = (RHO_A - RHO_W) * G                       # Winkler modulus Δρg [N/m³]
    EP = E / (1 - NU**2); D = EP * H**3 / 12.0
    BETA = (K / (4 * D))**0.25; ALPHA = 1.0 / BETA
    m = Model("data/suite1_strength/elastic_deep_60km_V4")

    # FE profiles on the (undeformed, total-Lagrangian) mesh — MATERIAL-frame S, see the docstring
    x = m.x.copy()                                 # distance from the loaded edge [m]
    w_fe = m.topography()                          # deflection [m], positive down
    Q_fe = m.resultant_material("sigma_xz [Pa]")   # ∫S_xz dz on the reference section [N/m]  (NOT m.V(), which is Cauchy)
    P = abs(Q_fe[np.argmin(np.abs(m.xkm - trench_ref_km(m)))])   # applied end shear inferred from FE Q at the trench

    # analytic Hetényi semi-infinite beam, end shear P at x'=0
    w_an = (2 * P * BETA / K) * np.exp(-BETA * x) * np.cos(BETA * x)
    Q_an = -P * np.exp(-BETA * x) * (np.cos(BETA * x) - np.sin(BETA * x))

    # panel (c): the ELASTIC SHEAR PARABOLA — isolated with the constant-shear idealized beam (σzz≈0, no
    # foundation), where the shear is a pure parabola.  Normalized (τ/τ_max vs ζ), so absolute scale is irrelevant.
    mb = Model("data/idealized_beam_elastic")
    Hb = mb.z.max(); zeta = (mb.z - Hb / 2) / (Hb / 2)
    tau_b = mb.column(mb.array("sigma_xz [Pa]"), 0.5 * mb.xkm.max())  # clean mid-span column; material-frame S (docstring)
    tau_bn = tau_b / tau_b[np.argmin(np.abs(zeta))]                    # normalize by the mid-depth (peak) value → +1 at centre
    par_n = 1 - zeta**2                                                # normalized parabola

    # --- numbers for SI Text S1 ---
    sel = x <= 4 * ALPHA
    relerr_w = np.nanmax(np.abs(w_fe[sel] - w_an[sel])) / np.nanmax(np.abs(w_an))
    relerr_Q = np.nanmax(np.abs(Q_fe[sel] - Q_an[sel])) / np.nanmax(np.abs(Q_an))
    off0 = (w_fe[0] / w_an[0] - 1) * 100
    relerr_tau = np.nanmax(np.abs(tau_bn - par_n))
    print(f"P (end shear, from FE)     = {P/1e12:.3f} TN/m")
    print(f"α = 1/β = {ALPHA/1e3:.1f} km   h/α = {H/ALPHA:.3f}")
    print(f"w(0):  FE {w_fe[0]:.0f} m   analytic(EB) {w_an[0]:.0f} m   offset {off0:+.1f}%  (thick-plate, expected)")
    print(f"max |w_FE − w_an| / max|w_an|  over 0–4α = {100*relerr_w:.1f}%")
    print(f"max |Q_FE − Q_an| / max|Q_an|  over 0–4α = {100*relerr_Q:.1f}%")
    print(f"shear parabola (idealized beam): max rel. dev. = {100*relerr_tau:.2f}%")

    # --- acceptance: the SI numbers, asserted (not just printed — the 2026-09-11 regression slipped past a print) ---
    P_APPLIED = 4.0e12                                                   # end shear of elastic_deep_60km_V4 [N/m]
    assert relerr_tau <= 0.005, f"shear parabola misfit {100*relerr_tau:.2f}% > 0.5% — wrong stress frame?"
    assert relerr_Q <= 0.01,    f"V misfit {100*relerr_Q:.2f}% > 1%"
    assert abs(P / P_APPLIED - 1) <= 0.01, f"end resultant {P/1e12:.3f} TN/m not within 1% of {P_APPLIED/1e12:.0f} TN/m"
    assert 0.005 <= off0 / 100 <= 0.03, f"thick-plate offset {off0:+.1f}% outside the expected 0.5–3% band"
    print("benchmark assertions passed")
    write_table("benchmark_boef", ["quantity", "value", "unit"], [
        ("end_resultant_from_FE", P / 1e12, "TN/m"), ("applied_end_shear", P_APPLIED / 1e12, "TN/m"),
        ("alpha", ALPHA / 1e3, "km"), ("h_over_alpha", H / ALPHA, "-"), ("w0_FE", float(w_fe[0]), "m"), ("w0_thin_beam", float(w_an[0]), "m"),
        ("w0_offset_vs_thin_beam", off0, "%"),
        ("w_misfit_max_0_4alpha", 100 * relerr_w, "%"), ("V_misfit_max_0_4alpha", 100 * relerr_Q, "%"),
        ("shear_parabola_misfit_max", 100 * relerr_tau, "%")],
        script="scripts/render_benchmark.py", figure="figures/benchmark_boef.png",
        models=["data/suite1_strength/elastic_deep_60km_V4", "data/idealized_beam_elastic"])

    # --- figure ---
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)
    xk = x / 1e3; xmax = 4 * ALPHA / 1e3

    ax[0].plot(xk, w_an, color="#b2182b", lw=3.0, alpha=0.45, label="analytic (Hetényi)")
    ax[0].plot(xk, w_fe, color="#111", lw=1.4, ls="--", label="FE (finite-strain)")
    ax[0].axhline(0, color="0.8", lw=0.7); ax[0].set_xlim(0, xmax); ax[0].invert_yaxis()
    ax[0].set_xlabel("distance from load  $x'$  [km]"); ax[0].set_ylabel("deflection  $w$  [m] (down)")
    ax[0].set_title("(a)  deflection", fontsize=10.5)
    ax[0].annotate(f"$w(0)$: FE {off0:+.1f}% vs EB\n(thick-plate offset, $h/\\alpha$={H/ALPHA:.2f})",
                   (0.03 * xmax, w_fe[0]), fontsize=8, va="bottom", color="0.3")
    ax[0].legend(fontsize=8.5, loc="lower right"); ax[0].grid(alpha=0.2)

    ax[1].plot(xk, Q_an/1e12, color="#b2182b", lw=3.0, alpha=0.45, label="analytic")
    ax[1].plot(xk, Q_fe/1e12, color="#111", lw=1.4, ls="--", label="FE")
    ax[1].axhline(0, color="0.8", lw=0.7); ax[1].set_xlim(0, xmax)
    ax[1].set_xlabel("distance from load  $x'$  [km]"); ax[1].set_ylabel(r"shear resultant  $V=\int\tau_{zx}\,dz$  [TN/m]")
    ax[1].set_title("(b)  shear resultant", fontsize=10.5)
    ax[1].legend(fontsize=8.5); ax[1].grid(alpha=0.2)

    ax[2].plot(par_n, zeta, color="#b2182b", lw=3.0, alpha=0.45, label=r"parabola $1-\zeta^2$")
    ax[2].plot(tau_bn, zeta, color="#111", lw=1.4, ls="--", label="FE (idealized beam)")
    ax[2].axvline(0, color="0.8", lw=0.7); ax[2].set_ylim(1, -1)
    ax[2].set_xlabel(r"normalized shear  $\sigma_{xz}/\sigma_{xz}^{\max}$")
    ax[2].set_ylabel(r"normalized depth  $\zeta=(z-h/2)/(h/2)$")
    ax[2].set_title("(c)  shear-stress parabola", fontsize=10.5)
    ax[2].annotate(f"max dev {100*relerr_tau:.1f}%", (0.05, 0.8), fontsize=8, color="0.3")
    ax[2].legend(fontsize=8.5, loc="upper right"); ax[2].grid(alpha=0.2)

    fig.suptitle("Methodology benchmark: finite-strain FE vs analytic beam-on-elastic-foundation", fontsize=11.5)
    OUT = "figures/benchmark_boef.png"
    fig.savefig(OUT, dpi=140); print("wrote", OUT)


if __name__ == "__main__":
    main()
