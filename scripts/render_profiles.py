"""render_profiles.py — depth-profile companion to the hero figure (elastic / Tresca / depth-dependent
von Mises, DD-VM), MASSLESS models, in the TRUE DEFORMED frame via the deformed-Cauchy pipeline.

Three rows (one rheology each), four DEPTH PROFILES per row (value vs deformed depth), at the SAME
flexure stations as the hero — max M · mid · max V.  Every field is the CAUCHY stress interpolated onto
the deformed vertical line (gpe_analysis.Model.deformed_line); massless ⇒ stress is zero outside the
plate (no hydrostatic caps).  Panels:
  (a) (σxx−σzz)(z)  [MPa]    the N_D integrand
  (b) σxz(z)        [MPa]
  (c) ρ̂ = τzx,x/g   [kg/m³]  equivalent density (FE Cauchy-shear gradient / g, branch B); depth centroid ●
                            marked.  Toggle PANEL3="grad" (or env PROFILES_PANEL3=grad) to plot raw τzx,x [kPa/m].
  (d) −Δσzz(z)      [MPa]    FILLED: area = the trench pull ΔGPE* (referenced to the TRENCH column, brought
                            inboard to the leftmost complete deformed column, at the same absolute level).

    python scripts/render_profiles.py
"""
import sys, os; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from gpe_analysis import Model, reference_lines, trench_pull, trench_ref_km, deformed_shear_gradient, write_table

MODELS = [("data/suite1_strength/elastic_deep_60km_V4",    "elastic"),
          ("data/suite1_strength/tresca_deep_150_60km_V4",  "elasto-plastic\n(Tresca)"),
          ("data/suite1_strength/dd_vm_asym_60km_V4",    "DD-VM asym\n(weak top)"),
          ("data/suite1_strength/dd_vm_sym_60km_V4",    "DD-VM sym\n(weak surf.)")]
OUT = "figures/profiles.png"
STATION_COLORS = {"moment_max": "#c0392b", "mid": "#8e44ad", "shear_max": "#16a085"}
STATION_LABELS = {"moment_max": "max $M$", "mid": "mid", "shear_max": "first isostatic ($x_I$)"}
ORDER = ["moment_max", "shear_max"]
QT_OFFSET_KM = 0.0     # (quasi-)trench column = leftmost complete deformed column + this offset;
                       # raise it to bring the trench profile inboard, clear of the loaded-edge layer
GRAV = 9.81
# THIRD COLUMN toggle: paper-consistent EQUIVALENT DENSITY ρ̂ = τzx,x/g [kg m⁻³] by default; set PANEL3
# to "grad" (or env PROFILES_PANEL3=grad) to plot the raw shear-stress gradient τzx,x [kPa m⁻¹] instead.
# Both are the same field up to the constant 1/g, so the centroid (●) is identical either way.
PANEL3 = os.environ.get("PROFILES_PANEL3", "rho")
P3SCALE, P3LABEL = {"rho":  (lambda t: t / GRAV, r"$\hat\rho=\tau_{zx,x}/g$  [kg m$^{-3}$]"),
                    "grad": (lambda t: t / 1e3,  r"$\tau_{zx,x}(z)$  [kPa m$^{-1}$]")}[PANEL3]


def stations(m):
    L = reference_lines(m)
    return {"moment_max": L["moment_max"],
            "mid": 0.5 * (L["moment_max"] + L["shear_max"]),
            "shear_max": L["shear_max"]}


DRHOG = (3300.0 - 1000.0) * 9.81      # water–rock contrast (only for the truncated above-plate GPE diagnostic)


def szz_on_grid(m, x_km, zgrid):
    """Cauchy σzz on the absolute-depth zgrid for the deformed column at x [km]; ZERO outside the plate.
    Massless ⇒ there genuinely is no stress above the deflected plate or below compensation, so the box
    integral needs NO caps — this is the force-balance integral itself.  The trench pull is this integral
    for the TRENCH column minus the ISOSTATIC column; the isostatic σzz≈0 (w=0 ⇒ no spring load), so the
    difference lives only where the plates are, and it is single-signed."""
    L = m.deformed_line(x_km); w, base = L["z"][0], L["z"][-1]
    return np.where((zgrid >= w) & (zgrid <= base), np.interp(zgrid, L["z"], L["szz"]), 0.0)


def plot_row(m, axes, is_top, zgrid):
    xst = stations(m)
    x_tr = trench_ref_km(m)
    szz_iso = szz_on_grid(m, xst["shear_max"], zgrid)                   # ISOSTATIC reference column (σzz ≈ 0)
    zt = zgrid / 1e3
    for k, (ax, (lab, kind)) in enumerate(zip(axes, [
            (r"$(\sigma_{xx}-\sigma_{zz})(z)$  [MPa]", "n_d"),
            (r"$\sigma_{xz}(z)$  [MPa]",               "sxz"),
            (P3LABEL,                                  "tau"),
            (r"$\sigma_{zz}(z)-\sigma_{zz}^{\rm iso}(z)$  [MPa]", "dszz")])):
        if kind == "dszz":
            # Single definition: σzz between fixed levels, zero outside plate, each column MINUS the isostatic
            # column. Trench−isostatic (the pull) and every intermediate column are single-signed. Area of the
            # TRENCH curve = the trench pull (= trench_pull ΔGPE*).  No caps.
            for plab, xk, col, lw, fill in [("trench", x_tr, "k", 2.6, True),
                                            ("max $M$", xst["moment_max"], STATION_COLORS["moment_max"], 1.6, False)]:
                prof = (szz_on_grid(m, xk, zgrid) - szz_iso) / 1e6
                if fill:
                    ax.fill_betweenx(zt, 0, prof, color=col, alpha=0.12)
                ax.plot(prof, zt, color=col, lw=lw, label=plab)
        else:
            # PLATE-SURFACE FRAME for panels (a)-(c): each curve is plotted against depth below THIS column's
            # own deflected top (z − z_top), so the rigid deflection is undone and the stress SHAPES align.
            for key in ORDER:
                xk = xst[key]; col = STATION_COLORS[key]
                if kind == "tau":
                    L0, tau = deformed_shear_gradient(m, xk)      # raw FE Cauchy-shear gradient (no smoothing)
                    ax.plot(P3SCALE(tau), (L0["z"] - L0["z"][0]) / 1e3, color=col, lw=2.0, label=STATION_LABELS[key])
                    tot = np.trapz(tau, L0["z"])                        # centroid of τzx,x (skip where ∫≈0, i.e. max V)
                    if np.abs(tot) > 0.2 * np.trapz(np.abs(tau), L0["z"]):
                        zc = np.trapz(L0["z"] * tau, L0["z"]) / tot
                        ax.plot(P3SCALE(np.interp(zc, L0["z"], tau)), (zc - L0["z"][0]) / 1e3,
                                marker="o", ms=8, color=col, mec="k", mew=0.7, zorder=6)
                else:
                    L0 = m.deformed_line(xk); prof = (L0["n_d"] if kind == "n_d" else L0["sxz"]) / 1e6
                    ax.plot(prof, (L0["z"] - L0["z"][0]) / 1e3, color=col, lw=2.0, label=STATION_LABELS[key])
            # the (quasi-)trench column (brought inboard by QT_OFFSET_KM); in the τzx,x panel its lobes +
            # centroid (●) are the equivalent-density dipole.  Panel (d) omits it (it is a difference vs the
            # isostatic column) and stays in the datum frame.
            x_qt = trench_ref_km(m) + QT_OFFSET_KM
            if kind == "tau":
                Lq, tauq = deformed_shear_gradient(m, x_qt)
                ax.plot(P3SCALE(tauq), (Lq["z"] - Lq["z"][0]) / 1e3, color="k", lw=2.6, label="trench", zorder=5)
                totq = np.trapz(tauq, Lq["z"])
                if np.abs(totq) > 0.2 * np.trapz(np.abs(tauq), Lq["z"]):
                    zcq = np.trapz(Lq["z"] * tauq, Lq["z"]) / totq
                    ax.plot(P3SCALE(np.interp(zcq, Lq["z"], tauq)), (zcq - Lq["z"][0]) / 1e3,
                            marker="o", ms=8, color="k", mec="w", mew=0.8, zorder=7)
            else:
                Lq = m.deformed_line(x_qt); profq = (Lq["n_d"] if kind == "n_d" else Lq["sxz"]) / 1e6
                ax.plot(profq, (Lq["z"] - Lq["z"][0]) / 1e3, color="k", lw=2.6, label="trench", zorder=5)
        ax.axvline(0, color="0.6", lw=0.7); ax.axhline(0, color="0.75", lw=0.6, ls=":")   # sea level
        ax.grid(alpha=0.22)
        (ax.set_title if is_top else ax.set_xlabel)(lab, fontsize=11)
    # dashed plate mid-depth reference: h/2 in the plate-surface frame (panels a-c); the absolute deformed
    # mid-depth in the datum frame (panel d).
    zmid = np.mean([0.5 * (m.deformed_line(xst[k])["z"][0] + m.deformed_line(xst[k])["z"][-1]) for k in ORDER])
    for j, ax in enumerate(axes):
        yline = zmid / 1e3 if j == 3 else (m.H / 2 / 1e3)
        ax.axhline(yline, color="k", lw=1.9, ls=(0, (6, 3)), alpha=0.8, zorder=1,
                   label=("plate mid-depth ($h/2$)" if (is_top and j == 0) else None))
    if is_top:
        axes[3].legend(fontsize=8.5, loc="lower right", title="column", framealpha=0.9)
    # annotate the trench pull (the trench-curve area) — the ONE definition, trench vs isostatic, no caps
    dGPE = trench_pull(m)[0]
    axes[3].text(0.055, 0.05, f"trench area $=$ pull\n$\\Delta$GPE$^* = {dGPE/1e12:+.2f}$ TN m$^{{-1}}$",
                 transform=axes[3].transAxes, fontsize=9.5, va="bottom", fontweight="bold", color="k",
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.9))


def main():
    ms = [Model(d) for d, _ in MODELS]
    ztop = min(float((m.z[0] + m.array("u", 1)[:, 0]).min()) for m in ms)     # shallowest deformed top (forebulge up)
    zbot = max(float((m.z[-1] + m.array("u", 1)[:, -1]).max()) for m in ms)   # deepest deformed base
    zgrid = np.linspace(ztop, zbot, 500)
    fig, axes = plt.subplots(len(MODELS), 4, figsize=(13.0, 3.55 * len(MODELS)), sharex="col", sharey=False, constrained_layout=True)
    print("panel-(d) self-check: trench-curve area (trench−isostatic, zero outside) == trench_pull ΔGPE*:")
    rows = []
    for i, (m, (_, rlab)) in enumerate(zip(ms, MODELS)):
        plot_row(m, axes[i], is_top=(i == 0), zgrid=zgrid)
        axes[i, 0].set_ylabel(f"{rlab}\n\ndepth below plate surface  [km]", fontsize=12)
        # self-check: the ONE definition — ∫(σzz_trench − σzz_iso) over the fixed box == trench_pull
        x_tr = trench_ref_km(m)
        area = np.trapz(szz_on_grid(m, x_tr, zgrid) - szz_on_grid(m, stations(m)["shear_max"], zgrid), zgrid)
        dGPE = trench_pull(m)[0]; wterm = DRHOG * m.deformed_line(x_tr)["z"][0] ** 2 / 2
        rows.append((MODELS[i][0], rlab.replace(chr(10), " "), area / 1e12, dGPE / 1e12, 100 * abs(area - dGPE) / abs(dGPE)))
        print(f"   {rlab.replace(chr(10),' '):24s}: area={area/1e12:+.4f}  trench_pull={dGPE/1e12:+.4f}  (Δ={100*abs(area-dGPE)/abs(dGPE):.2f}%)   [truncated above-plate GPE ≈ {wterm/1e12:.3f}]")
    Hkm = ms[0].H / 1e3                                                      # (a)-(c): plate-surface frame (0 = plate top)
    for ax in axes[:, :3].flat:
        ax.set_ylim(Hkm + 2, -2)
    for ax in axes[:, 3]:                                                     # (d): datum frame (absolute z), axis on the right
        ax.set_ylim(zbot / 1e3, ztop / 1e3)
        ax.yaxis.set_label_position("right"); ax.yaxis.tick_right()
    for ax in list(axes[:, 1]) + list(axes[:, 2]):                           # hide interior duplicate y-ticklabels
        ax.tick_params(labelleft=False)
    axes[0, 3].set_ylabel("depth below reference level  [km]", fontsize=10)
    axes[0, 0].legend(fontsize=9.5, loc="lower left", frameon=True, framealpha=0.9, title="location")
    fig.suptitle("Depth profiles at the flexure locations (trench · max $M$ · first isostatic $x_I$); (a)-(c) below plate surface, (d) below reference level   "
                 "● = centroid of $\\tau_{zx,x}$", fontsize=12.5)
    fig.savefig(OUT, dpi=150); print("wrote", OUT)
    write_table("profiles_selfcheck", ["model", "label", "trench_curve_area_TN", "trench_pull_TN", "diff_pct"], rows,
                script="scripts/render_profiles.py", figure="figures/profiles.png", models=[d for d, _ in MODELS])


if __name__ == "__main__":
    if len(sys.argv) > 1:        # optional: render to a preview path instead of the committed figure
        OUT = sys.argv[1]
    main()
