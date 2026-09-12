"""render_thickness_compare.py — SI Suite 4: does plate THICKNESS control the trench pull?

Four uniform-Tresca models (σ_Y = 150 MPa, identical strength profile), each tuned to the SAME trench
deflection (3233 m) but a different mechanical thickness h ∈ {30, 40, 50, 60} km — h=60 is the locked
baseline, {30,40,50} are the suite4_thickness sweep (V bisected per h; see paper_models.jl `thickness`).
Only h changes; the first isostatic column is located separately for each (it moves inboard with the
flexural wavelength).  Message: with the finite thickness set, the strength-profile details barely matter
(Fig 5), but the THICKNESS itself enters — a thinner plate supports the topography over a shallower dipole
(arm ~h/2), so the same trench deflection yields a smaller pull, roughly pull ∝ h.

  (a) pull ΔGPE* vs applied load V (bottom) with thickness h on a twin TOP axis (V∝h) — one line carries
      pull ∝ load ∝ thickness; through-origin fit ΔGPE* ≈ 0.65 V.  (Merges the former thickness_paths panel.)
  (b) deflection profiles, all matched at the trench (referenced to 0 at each model's isostatic column).
  (c) equivalent density ρ̂ = τzx,x/g at the trench (plate-surface frame; ● centroid = dipole arm; h/2 ref).
  (d) equivalent density at the max-moment (hinge) column.

  /opt/anaconda3/envs/pyvista-env/bin/python render_thickness_compare.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from gpe_analysis import Model, trench_ref_km, reference_lines, trench_pull, deformed_shear_gradient

#           model dir                                             h[km]  colour (thin→warm, thick→cool)
MODELS = [("data/suite4_thickness/tresca_150_30km",             30,   "#c0392b"),
          ("data/suite4_thickness/tresca_150_40km",             40,   "#e67e22"),
          ("data/suite4_thickness/tresca_150_50km",             50,   "#2779b5"),
          ("data/suite1_strength/tresca_deep_150_60km_V4",      60,   "#1f3b73")]
VLOAD = {30: 2.045, 40: 2.702, 50: 3.354, 60: 4.000}   # tuned load per plate [TN/m] (tuned_V.txt; 60 = baseline)
OUT = "data/reference_figures/thickness_compare.png"
WINDOW_KM = 220
GRAV = 9.81


def load(mdir):
    m = Model(mdir)
    R = reference_lines(m)                                                  # shear_max = first isostatic column; moment_max = hinge
    x_tr = trench_ref_km(m)
    xs = np.linspace(x_tr, WINDOW_KM, 90)
    Szz0 = m.deformed_resultants(x_tr)[0]
    dgpe = np.array([-(m.deformed_resultants(x)[0] - Szz0) for x in xs]) / 1e12
    return dict(m=m, x_tr=x_tr, xi=R["shear_max"], x_mm=R["moment_max"], xs=xs,
                dgpe=dgpe, pull=trench_pull(m)[0] / 1e12, H=m.H)


def _ylim(ax, data):
    ax.set_ylim(max(d["H"] for d in data) / 1e3 + 2, -2); ax.axvline(0, color="0.6", lw=.7)


def _rho(ax, data, xkey, title):
    """Equivalent density ρ̂ = τzx,x/g at column `xkey`, plate-surface frame; centroid ● + h/2 per plate."""
    for (mdir, hk, col), d in zip(MODELS, data):
        m = d["m"]; L, tau = deformed_shear_gradient(m, d[xkey]); z = (L["z"] - L["z"][0]) / 1e3
        ax.plot(gaussian_filter1d(tau / GRAV, 1.5), z, color=col, lw=2.2, label=f"$h={hk}$ km")   # smooth yield-front spikes
        # signed charge centroid z_c = ∫z·τ dz / ∫τ dz (the dipole depth) from RAW τ.  Divides by the net
        # charge, so it is only meaningful when the profile is charge-DOMINATED; when |∫τ| is small vs ∫|τ|
        # (a near-balanced dipole, e.g. the thin-plate hinge) the centroid is ill-conditioned — suppress it.
        tot = np.trapz(tau, L["z"]); absq = np.trapz(np.abs(tau), L["z"])
        if abs(tot) > 0.45 * absq:
            zc = (np.trapz(L["z"] * tau, L["z"]) / tot - L["z"][0]) / 1e3
            ax.plot(0.0, zc, "o", color=col, ms=9, mec="k", mew=0.9, zorder=6)   # centroid DEPTH, marked on the zero line (not on the curve)
        ax.axhline(d["H"] / 2 / 1e3, color=col, lw=1.0, alpha=.4)           # h/2 reference (thin solid)
    _ylim(ax, data)
    ax.set_xlabel(r"$\hat\rho = \tau_{zx,x}/g$  [kg m$^{-3}$]"); ax.set_ylabel("depth below plate surface  [km]")
    ax.set_title(title, fontsize=11); ax.legend(fontsize=8.5, ncol=2); ax.grid(alpha=.25)


def main():
    data = [load(d) for d, _, _ in MODELS]
    hs = np.array([hk for _, hk, _ in MODELS], float)
    pulls = np.array([d["pull"] for d in data])
    fig, axf = plt.subplots(2, 2, figsize=(12.5, 10.6), constrained_layout=True); ax = axf.flat

    # (a) CONDENSED: pull vs applied load V (bottom axis) with plate thickness h on a twin TOP axis — the two
    # quantities that vary across the suite.  The top ticks sit on the same four points because V ∝ h (V/h
    # constant to ~2%), so one line carries pull ∝ load ∝ thickness; the through-origin fit is ΔGPE* ≈ 0.65 V.
    Vs = np.array([VLOAD[hk] for _, hk, _ in MODELS])
    s0 = float(np.sum(Vs * pulls) / np.sum(Vs**2))                          # proportional (through-origin) slope
    xl = (1.35, Vs.max() * 1.06)                                            # zoom onto the data (h≥20 km); origin not shown
    ax[0].plot([xl[0], xl[1]], [s0 * xl[0], s0 * xl[1]], "-", color="0.55", lw=1.4, zorder=1,
               label=rf"$\Delta$GPE$^*\approx{s0:.2f}\,V$  (through origin)")
    for (mdir, hk, col), d in zip(MODELS, data):
        ax[0].plot(VLOAD[hk], d["pull"], "o", color=col, ms=12, mec="k", mew=1.0, zorder=5, label=f"$h={hk}$ km")
    ax[0].set_xlim(*xl); ax[0].set_ylim(pulls.min() - 0.28, pulls.max() + 0.18)
    ax[0].set_xlabel(r"applied load  $V$  (end shear force)  [TN m$^{-1}$]")
    ax[0].set_ylabel(r"trench pull  $\Delta$GPE$^{*}$  [TN m$^{-1}$]")
    ax[0].set_title("(a)  pull vs load and thickness  (matched deflection)", fontsize=11)
    ax[0].legend(fontsize=8.5, loc="upper left"); ax[0].grid(alpha=.25)
    axT = ax[0].twiny(); axT.set_xlim(*xl)                                  # thickness on top: ticks at each model (V∝h)
    axT.set_xticks(Vs); axT.set_xticklabels([f"{hk}" for _, hk, _ in MODELS])
    axT.set_xlabel(r"plate thickness  $h$  [km]  ($V\propto h$)", fontsize=9.5)

    # (b) deflection vs distance (referenced to zero at each model's first isostatic column) — the matched input
    for (mdir, hk, col), d in zip(MODELS, data):
        xk = d["m"].xkm - d["x_tr"]; w = d["m"].topography() / 1e3
        w = w - np.interp(d["xi"] - d["x_tr"], xk, w)
        sel = (xk >= 0) & (xk <= WINDOW_KM - d["x_tr"])
        ax[1].plot(xk[sel], w[sel], color=col, lw=2.2, label=f"$h={hk}$ km")
        ax[1].plot(d["xi"] - d["x_tr"], 0.0, "o", color=col, ms=7, mec="k", mew=0.8, zorder=6)   # isostatic column (w=0)
    ax[1].invert_yaxis(); ax[1].axhline(0, color="0.7", lw=0.7)            # deflection positive-down ⇒ trench dips
    ax[1].set_xlabel("distance from trench  [km]"); ax[1].set_ylabel(r"deflection  $w$  [km]  ($=0$ at isostatic)")
    ax[1].set_title("(b)  deflection  (matched at the trench; ● first isostatic column)", fontsize=11)
    ax[1].legend(fontsize=8.5, loc="lower right", ncol=2); ax[1].grid(alpha=.25)

    # (c,d) equivalent density — the dipole (thinner plate: higher amplitude, shallower)
    _rho(ax[2], data, "x_tr", "(c)  equivalent density at the trench")
    _rho(ax[3], data, "x_mm", "(d)  equivalent density at max moment")
    # h=30 max-moment centroid: the hinge is fully plastic and the signed shear-gradient centroid is ill-defined
    # (net dipole ≈ 0), so it is SUPPRESSED above.  Plot an OPEN circle at the depth INTERPOLATED from the three
    # well-resolved thicker plates (a linear centroid-depth vs h fit) — NOT from the h=30 gradient.  Caption flags it.
    mm_h, mm_z = [], []
    for (mdir, hk, col), d in zip(MODELS, data):
        if hk == 30:
            continue
        L, tau = deformed_shear_gradient(d["m"], d["x_mm"]); z = L["z"]
        mm_h.append(hk); mm_z.append((np.trapz(z * tau, z) / np.trapz(tau, z) - z[0]) / 1e3)
    z30 = np.polyval(np.polyfit(mm_h, mm_z, 1), 30)
    ax[3].plot(0.0, z30, "o", mfc="none", mec=MODELS[0][2], mew=2.0, ms=11, zorder=7,
               label=f"$h=30$ (interp., {z30:.0f} km)")
    ax[3].legend(fontsize=8.5, ncol=2)
    print(f"  h=30 max-moment centroid (interpolated from h=40/50/60): {z30:.1f} km")

    fig.suptitle(r"Plate thickness controls the trench pull  (uniform Tresca, matched deflection, $\sigma_Y=150$ MPa)", fontsize=13.5)
    outp = sys.argv[1] if len(sys.argv) > 1 else OUT
    fig.savefig(outp, dpi=150); print("wrote", outp)
    print(f"  through-origin ΔGPE*/V = {s0:.3f};  ratios = " + ", ".join(f"{d['pull']/VLOAD[hk]:.3f}" for (mdir, hk, col), d in zip(MODELS, data)))
    for (mdir, hk, col), d in zip(MODELS, data):
        w = d["m"].deformed_line(d["x_tr"])["z"][0]
        print(f"  h={hk:2d} km: w={w/1e3:.2f} km  iso {d['xi']-d['x_tr']:.0f} km  maxM {d['x_mm']-d['x_tr']:.0f} km  pull {d['pull']:.2f}  h/2={d['H']/2/1e3:.0f}")


if __name__ == "__main__":
    main()
