"""render_corrected_density.py — EXPLORATORY: the trench pull as a dipole of the CORRECTED density.

Corrected density  ρ̂_c = true density (water above the deflected top · plate 3300 · mantle 3300) + equivalent
density (shear support τ_zx,x/g in the plate).  Three columns — trench · max M · isostatic (the reference).

  (1) corrected density  ρ̂_c(z)          — symlog(ρ−3300): the huge water↔rock step AND the mild equivalent-density wiggle
  (2) cumulative average ρ̄(z)=z⁻¹∫ρ̂_c    — every column → the SAME constant (equal integrated mass = isostasy)
  (3) the dipole                          — a-priori split: TRUE (water @ w/2) ↔ equivalent (shear); g·m·arm ≈ pull

Note on the arm (~34 km, > h/2): the equivalent-density lobe sits below the deflected mid-plate because (a) the deformed
line's top is the deflected depth (~+3 km) and (b) the trench-zone shear support is genuinely bottom-heavy
(plate-relative centroid ~33 km, steady over x=5-40 km).  This ~34 km arm is what reproduces the measured pull
(g·Δρ·w·arm ≈ 2.5 TN/m) and matches the manuscript scaling ½(w_T+z_np); h/2 would under-predict by ~15%.

    python scripts/render_corrected_density.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.transforms import blended_transform_factory
from gpe_analysis import Model, reference_lines, trench_ref_km, deformed_shear_gradient, trench_pull, write_table, RHO_M, RHO_W, G, CENTROID_RATIO_MIN

m = Model("data/suite1_strength/tresca_deep_150_60km_V4")
RW, RP = RHO_W, RHO_M; DR = RP - RW
rl = reference_lines(m); xt = trench_ref_km(m); xmm = rl["moment_max"]; xi = rl["isostatic"]
COLS = [("trench", xt, "#b2182b"), (r"max $M$", xmm, "#ef8a62"), ("isostatic", xi, "#111111")]
ZG = np.linspace(0, 90e3, 9001); DZ = ZG[1] - ZG[0]

_Lt = deformed_shear_gradient(m, xt)[0]["z"]                      # trench deformed column geometry
W_T, BASE_T = _Lt[0], _Lt[-1]; MID_T = 0.5 * (W_T + BASE_T)      # deflected top / base / mid-plate
YTOP, YBOT = -3.0, BASE_T / 1e3 + 12.0                           # shared depth window — a little uniform mantle below


def reflines(a):                                                 # deflected mid-plate & base, shared across panels
    a.axhspan(BASE_T / 1e3, YBOT, color="0.90", zorder=0)        # uniform mantle below the deflected base
    a.axhline(MID_T / 1e3, color="0.55", lw=0.9, ls=(0, (6, 4)), zorder=1)
    a.axhline(BASE_T / 1e3, color="0.30", lw=1.2, zorder=1)
    a.set_ylim(YBOT, YTOP)


def corrected(xk):                                               # full corrected-density profile on ZG
    L, tau = deformed_shear_gradient(m, xk); z = L["z"]; top, base = z[0], z[-1]
    r = np.full_like(ZG, RP); r[ZG < top] = RW                   # water above the deflected top, else 3300
    inpl = (ZG >= top) & (ZG <= base); r[inpl] = RP + np.interp(ZG[inpl], z, tau / G)
    return r


def dipole_lobes(xk):                                            # a-priori: TRUE (water @ w/2) + equivalent (shear centroid)
    L, tau = deformed_shear_gradient(m, xk); z = L["z"]; w = z[0]
    m_true = DR * w; z_true = w / 2.0                            # water deficit — KNOWN from the deflected depth
    net = np.trapz(tau, z); m_ps = net / G                                          # shear support — from the FE gradient
    # signed centre of mass; undefined where the net charge vanishes (the isostatic column: |∫τ| ≪ ∫|τ|) → NaN, as every
    # other caller guards it (thickness_compare, hero, profiles, the notebook)
    z_ps = np.trapz(z * tau, z) / net if abs(net) > CENTROID_RATIO_MIN * np.trapz(np.abs(tau), z) else float("nan")
    return dict(w=w, m_true=m_true, z_true=z_true, m_ps=m_ps, z_ps=z_ps)


def main():
    rc = {lab: corrected(xk) for lab, xk, _ in COLS}
    zk = ZG / 1e3
    fig, ax = plt.subplots(1, 3, figsize=(12.4, 6.2), gridspec_kw={"width_ratios": [1, 1, 0.62]},
                           constrained_layout=True)

    # (1) corrected density, symlog about the 3300 baseline so water step AND equivalent-density wiggle both read
    for lab, xk, c in COLS:
        ax[0].plot(rc[lab] - RP, zk, color=c, lw=2.0, ls=("--" if lab == "isostatic" else "-"), label=lab)
    ax[0].set_xscale("symlog", linthresh=400.0)
    ticks = [-2300, -400, 0, 400]                               # 1000(water) · 2900 · 3300(rock) · 3700
    ax[0].xaxis.set_major_locator(FixedLocator(ticks))
    ax[0].xaxis.set_major_formatter(FixedFormatter([f"{RP + t:.0f}" for t in ticks]))
    ax[0].axvline(0, color="0.7", lw=0.8, zorder=0)
    ax[0].set_xlim(-2600, 700); reflines(ax[0])
    ax[0].set_xlabel(r"corrected density  $\hat\rho_c$  [kg m$^{-3}$]  (symlog about 3300)")
    ax[0].set_ylabel("depth  [km]")
    ax[0].set_title(r"(a)  corrected density  $\hat\rho_c=\rho+g^{-1}\partial_x\sigma_{xz}$", fontsize=10.5)
    ax[0].annotate("water deficit", (-900, -0.2), fontsize=8, ha="center", va="bottom", color="#b2182b")
    ax[0].annotate("equivalent\n(shear support)", (172, 55), fontsize=8, ha="left", va="center", color="0.3")
    tf = blended_transform_factory(ax[0].transAxes, ax[0].transData)
    ax[0].text(0.02, MID_T / 1e3, "deflected mid-plate", transform=tf, ha="left", va="bottom", fontsize=7.5, color="0.45")
    ax[0].text(0.02, BASE_T / 1e3, "deflected base", transform=tf, ha="left", va="bottom", fontsize=7.5, color="0.25")
    ax[0].legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 0.82)); ax[0].grid(alpha=0.2)

    # (2) cumulative average corrected density -> a constant with depth (equal integrated mass)
    for lab, xk, c in COLS:
        avg = np.cumsum(rc[lab]) * DZ / np.maximum(ZG, DZ)
        ax[1].plot(avg, zk, color=c, lw=2.0, ls=("--" if lab == "isostatic" else "-"), label=lab)
    conv = np.mean([np.cumsum(rc[l])[-1] * DZ / ZG[-1] for l, _, _ in COLS])
    ax[1].axvline(conv, color="0.5", lw=0.8, ls=":")
    ax[1].set_xlim(900, 3450); reflines(ax[1])
    ax[1].set_xlabel(r"cumulative average  $\bar\rho(z)=z^{-1}\!\int_0^z\!\hat\rho_c\,dz'$  [kg m$^{-3}$]")
    ax[1].set_ylabel("depth  [km]")
    ax[1].set_title("(b)  every column $\\to$ the same integrated density", fontsize=10.5)
    ax[1].annotate(f"$\\bar\\rho\\to${conv:.0f}\n(equal mass)", (conv - 120, 58), fontsize=8.5, ha="right",
                   va="bottom", style="italic", color="0.3")
    ax[1].grid(alpha=0.2)                                         # no legend — same columns as panel (1)

    # (3) the dipole = ΔGPE of each column RELATIVE TO THE ISOSTATIC COLUMN.  TRUE (water @ w/2) <-> equivalent (shear).
    Szz_i = m.deformed_resultants(xi)[0]                          # isostatic reference σzz-integral
    dgpe = lambda xk: -(Szz_i - m.deformed_resultants(xk)[0]) / 1e12     # ΔGPE vs isostatic [TN/m]
    glyph = COLS[:2]; slot = {lab: 1.0 + i for i, (lab, _, _) in enumerate(glyph)}   # trench · max M
    m0 = dipole_lobes(xt)["m_ps"]
    for lab, xk, c in glyph:
        d = dipole_lobes(xk); x = slot[lab]
        st = 2500 * abs(d["m_true"]) / m0; sp = 2500 * abs(d["m_ps"]) / m0
        arm = (d["z_ps"] - d["z_true"]) / 1e3
        ax[2].plot([x, x], [d["z_true"] / 1e3, d["z_ps"] / 1e3], color=c, lw=2.5, zorder=2)
        ax[2].scatter([x], [d["z_true"] / 1e3], s=st, c="#4292c6", ec="k", lw=0.6, zorder=3)   # TRUE (water) @ w/2
        ax[2].scatter([x], [d["z_ps"] / 1e3], s=sp, c=c, ec="k", lw=0.6, zorder=3)             # equivalent (shear), column colour
        off = (sp / 3.1416) ** 0.5 + 6
        ax[2].annotate(f"$\\Delta$GPE = {dgpe(xk):.2f}\nTN/m\n(arm {arm:.0f} km)", (x, d["z_ps"] / 1e3),
                       xytext=(off, 0), textcoords="offset points", va="center", fontsize=8, color=c)
    # identify the two lobe types, centred in black on the (large) trench dipole circles — replaces a legend
    dt = dipole_lobes(xt)
    ax[2].annotate("water\n(true)", (slot["trench"], dt["z_true"] / 1e3), ha="center", va="center",
                   fontsize=7.5, color="k", zorder=4)
    ax[2].annotate("shear\n(equivalent)", (slot["trench"], dt["z_ps"] / 1e3), ha="center", va="center",
                   fontsize=7.5, color="k", zorder=4)
    ax[2].axhline(0, color="0.85", lw=0.7)
    reflines(ax[2]); ax[2].set_xlim(0.45, 2.9)
    ax[2].set_xticks(list(slot.values())); ax[2].set_xticklabels([l for l, _, _ in glyph], fontsize=10)
    ax[2].set_ylabel("depth  [km]"); ax[2].set_xlabel(r"column   ($\Delta$GPE relative to isostatic)")
    ax[2].set_title(r"(c)  $\Delta$GPE dipole (vs isostatic column)", fontsize=10.5)
    ax[2].grid(alpha=0.2, axis="y")

    dg, _, _ = trench_pull(m)
    fig.suptitle(r"Corrected density: same integrated mass in every column $\Rightarrow$ trench pull is a density dipole"
                 f"  (measured pull {dg/1e12:.2f} TN/m)", fontsize=11.5)
    outp = sys.argv[1] if len(sys.argv) > 1 else "figures/corrected_density.png"
    rows = []
    for lab, xk, c in COLS:
        d = dipole_lobes(xk)
        rows.append((lab.replace("$", "").replace("\\", ""), xk, dgpe(xk), d["m_true"] * G / 1e6, d["z_true"] / 1e3, d["m_ps"] * G / 1e6, d["z_ps"] / 1e3, (d["z_ps"] - d["z_true"]) / 1e3))
    write_table("corrected_density", ["column", "x_km", "dGPE_vs_isostatic_TN", "water_lobe_pressure_MPa", "water_lobe_depth_km",
                                      "equivalent_lobe_pressure_MPa", "equivalent_lobe_depth_km", "dipole_arm_km"], rows,
                script="scripts/render_corrected_density.py", figure=outp, models=[m.dir],
                meta={"deflected_mid_plate_km": MID_T / 1e3, "deflected_base_km": BASE_T / 1e3, "trench_depth_km": W_T / 1e3,
                      "cumulative_average_density_at_base_kg_m3": float(conv)})
    fig.savefig(outp, dpi=140); print("wrote", outp)


if __name__ == "__main__":
    main()
