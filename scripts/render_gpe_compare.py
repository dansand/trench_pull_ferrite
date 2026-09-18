"""render_gpe_compare.py — Figure 5: the trench pull IS the equivalent-density dipole (Suite 1, 2×2).

ΔGPE*(x) = −Δσ̄zz, the box integral of the Cauchy σzz between fixed levels (zero outside the plate),
trench-referenced — the SINGLE definition used everywhere (Model.deformed_resultants).  From vertical
equilibrium ∂σzz/∂z = −τzx,x (massless: no body force), σzz is the load of the EQUIVALENT DENSITY
ρ̂ = τzx,x/g (purely the horizontal gradient of the vertical shear stress — no real density).

THE ARM MATTERS.  Integrating −∫z·τ by parts gives −∫σzz (= ΔGPE*) PLUS a boundary term σzz(top)·w.
Measured from SEA LEVEL (z) that term is left in → a drift growing as w decays off the trench; measured
from the PLATE TOP (z−w) it vanishes → the dipole equals ΔGPE* (~2%).  The reconstruction drawn here uses
the plate-top arm; the sea-level error is printed (not drawn) so the SI-quoted numbers stay reproducible.

ONE figure, one panel per Suite-1 strength model (elastic | Tresca | DD-VM asym | DD-VM sym), three curves:
  ΔGPE*(x) direct (grey band) · ρ̂-dipole reconstruction, plate-top arm (black, traces the band) ·
  mid-plate approximation Δρg·(h/2)·(w_T−w) (black dashed; its error is in each panel title).

    python scripts/render_gpe_compare.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from gpe_analysis import Model, trench_ref_km, deformed_shear_gradient, write_table, DRHOG

S1 = "data/suite1_strength"
# The Tresca baseline is the reference model; SUITE1 is the full 2×2 rheology panel.  (DD-VM entries are
# whatever currently sits in suite1_strength — refreshed when the friction DD-VMs are promoted.)
SUITE1 = [(f"{S1}/elastic_deep_60km_V4",    "elastic"),
          (f"{S1}/tresca_deep_150_60km_V4", "Tresca (uniform)"),
          (f"{S1}/dd_vm_asym_60km_V4",   "DD-VM (asymmetric)"),
          (f"{S1}/dd_vm_sym_60km_V4",   "DD-VM (symmetric)")]
OUT_SUITE1 = "figures/gpe_compare_suite1.png"     # Figure 5: Suite-1 robustness (2×2) — the ONLY output
WINDOW_KM = 200
# B&W-robust styling (2026-09-11, user): distinguish by weight/style, not hue, so the figure reads
# identically when printed in greyscale.  The truth is a thick pale-GREY band; the reconstruction is a
# thin BLACK line that traces it (agreement, both continuous — no discrete markers); the mid-plate
# approx is BLACK DASHED (the outlier).
C_BAND = "0.68"         # ΔGPE* (direct ∫σzz): thick pale-grey band
C_LINE = "#111111"      # ρ̂ dipole reconstruction: thin black solid, riding on the band


def profile(model_dir):
    """ΔGPE*(x); the equivalent-density dipole with the arm from SEA LEVEL and from the PLATE TOP; and the
    mid-plate approximation — all trench-referenced, one definition (Cauchy σzz between fixed levels)."""
    m = Model(model_dir)
    x_ref = trench_ref_km(m) + 5.0; H = m.H                          # reference column: first with CENTRED τ (clean)
    xs = np.linspace(x_ref, WINDOW_KM, 58)
    Szz_r, _ = m.deformed_resultants(x_ref)
    Lr, tr = deformed_shear_gradient(m, x_ref); w_ref = Lr["z"][0]
    pdsea_r = -np.trapz(Lr["z"] * tr, Lr["z"])                       # arm from sea level (z=0)
    pdtop_r = -np.trapz((Lr["z"] - w_ref) * tr, Lr["z"])            # arm from plate top (z−w)
    dgpe, pdsea, pdtop, approx = [], [], [], []
    for xk in xs:
        Szz, _ = m.deformed_resultants(xk); dgpe.append(-(Szz - Szz_r))
        L, t = deformed_shear_gradient(m, xk); z = L["z"]           # centred everywhere on this window
        pdsea.append(-np.trapz(z * t, z) - pdsea_r)
        pdtop.append(-np.trapz((z - z[0]) * t, z) - pdtop_r)
        approx.append(DRHOG * (H / 2) * (w_ref - z[0]))             # mid-plate approx  Δρg·(h/2)·(w_ref−w)
    pdsea, pdtop = gaussian_filter1d(np.array(pdsea), 1.0), gaussian_filter1d(np.array(pdtop), 1.0)
    return m, xs, np.array(dgpe), pdsea, pdtop, np.array(approx)


def _plot_pull(ax, x, dgpe, pdtop, approx, s=1e12, legend=False):
    """Top-panel content: true ΔGPE*, the equivalent-density dipole (plate-top arm), the mid-plate approx."""
    ax.plot(x, dgpe / s, color=C_BAND, lw=6.0, solid_capstyle="round", zorder=2,
            label=r"$\Delta\mathrm{GPE}^{*}$ (direct $-\!\int\sigma_{zz}\,dz$)")
    ax.plot(x, pdtop / s, color=C_LINE, lw=1.5, solid_capstyle="round", zorder=4,
            label=r"$\hat\rho$ dipole reconstruction (plate-top arm)")
    ax.plot(x, approx / s, color=C_LINE, lw=1.5, ls=(0, (5, 3)), zorder=3,
            label=r"mid-plate approx  $\Delta\rho g\,\frac{h}{2}(w_T-w)$")
    ax.axhline(0, color="0.7", lw=0.6); ax.grid(alpha=0.2)
    if legend:
        ax.legend(fontsize=10.5, loc="upper left", handlelength=2.8, framealpha=0.9)


def _errs(dgpe, pdsea, pdtop, approx):
    d = np.max(np.abs(dgpe))
    return (100 * np.nanmax(np.abs(pdsea - dgpe)) / d, 100 * np.nanmax(np.abs(pdtop - dgpe)) / d,
            100 * np.max(np.abs(approx - dgpe)) / d)


def main_suite1():
    """Full Suite-1 robustness: pull = equivalent-density dipole across all four rheologies (2×2 top-panels)."""
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.4), sharex=True, sharey=True, constrained_layout=True)
    rows = []
    for k, (mdir, title) in enumerate(SUITE1):
        ax = axes.flat[k]
        m, x, dgpe, pdsea, pdtop, approx = profile(mdir)
        _plot_pull(ax, x, dgpe, pdtop, approx, legend=(k == 0))
        esea, etop, eap = _errs(dgpe, pdsea, pdtop, approx)
        # The plate-top-arm reconstruction error (etop, ~2%) and the sea-level-arm error (esea, ~8%) are the
        # numbers quoted in the SI; they are printed (not drawn) so the quoted values stay reproducible.
        # The h/2 miss in the title is the physical, reportable one.
        print(f"   {title:22s}: plate-top arm {etop:.1f}%,  sea-level arm {esea:.1f}%,  mid-plate approx {eap:.1f}%")
        rows.append((mdir, title, etop, esea, eap))
        ax.set_title(f"{title}    ($h/2$ estimate off {eap:.0f}%)", fontsize=10.5)
        if k % 2 == 0:
            ax.set_ylabel(r"$\Delta\mathrm{GPE}^{*}$  [TN m$^{-1}$]", fontsize=11)
        if k >= 2:
            ax.set_xlabel(r"distance from trench  [km]", fontsize=11)
    fig.suptitle(r"Trench pull $\Delta\mathrm{GPE}^{*}$ = equivalent-density dipole across Suite 1 — "
                 "relatively insensitive to the strength model (plate-top arm)", fontsize=12.5)
    fig.savefig(OUT_SUITE1, dpi=135); print("wrote", OUT_SUITE1)
    write_table("gpe_compare_reconstruction", ["model", "label", "plate_top_arm_err_pct", "sea_level_arm_err_pct", "mid_plate_approx_err_pct"],
                rows, script="scripts/render_gpe_compare.py", figure="figures/gpe_compare_suite1.png", models=[d for d, _ in SUITE1])


if __name__ == "__main__":
    main_suite1()
