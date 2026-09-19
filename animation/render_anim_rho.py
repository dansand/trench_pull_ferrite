"""animation/render_anim_rho.py — companion to render_anim.py (STANDALONE, DISPOSABLE).

Same 24 load-step frames (animation/frames/step_NN), no re-solving.  Shows the SHEAR-STRESS GRADIENT /
equivalent density rho_hat = tau_zx,x / g through the loading cycle:
  (left)  rho_hat as a FIELD on the deformed plate.  FIXED symmetric-log colour scale (linear near zero so
          the interior structure shows, compressed at the extremes so the strong trench column does not
          saturate) — fixed across frames so the amplitude visibly GROWS with load.
  (right) rho_hat(z) DEPTH PROFILES at the TRENCH and the MAX-BENDING-MOMENT columns — the equivalent-
          density dipoles forming, each with its centroid (dipole arm) marked, against h/2.

Reuses ../gpe_analysis.py READ-ONLY.  Delete animation/ to remove everything.

  python animation/render_anim_rho.py
"""
import sys, os, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import SymLogNorm
from scipy.ndimage import zoom
from PIL import Image
from gpe_analysis import Model, trench_pull, trench_ref_km, reference_lines, deformed_shear_gradient

HERE = os.path.dirname(__file__)
FRAMES = sorted(glob.glob(os.path.join(HERE, "frames", "step_*")))
PNGDIR = os.path.join(HERE, "frames", "_png_rho"); os.makedirs(PNGDIR, exist_ok=True)
OUT_GIF = os.path.join(HERE, "trench_pull_rhohat.gif")
OUT_MP4 = os.path.join(HERE, "trench_pull_rhohat.mp4")
WINDOW_KM = 260
WARP = 6.0
GRAV = 9.81
LINTHRESH = 20.0                                  # |rho_hat| below this is shown linearly [kg m^-3]
COLS = [("trench",  trench_ref_km, "k",       2.6),
        ("max $M$", lambda m: reference_lines(m)["moment_max"], "#c0392b", 2.0)]


def meta_V(d):
    try:
        return float(open(os.path.join(d, "meta.txt")).read().split("V_TN=")[1].split()[0])
    except Exception:
        return np.nan


def rho_profile(m, xk):
    """rho_hat(z) [kg m^-3] and depth-below-top [km] at column xk, plus its signed centroid depth (or nan)."""
    L, tau = deformed_shear_gradient(m, xk); rho = tau / GRAV
    z = (L["z"] - L["z"][0]) / 1e3
    tot = np.trapz(tau, L["z"])
    zc = (np.trapz(L["z"] * tau, L["z"]) / tot - L["z"][0]) / 1e3 if abs(tot) > 0.2 * np.trapz(np.abs(tau), L["z"]) else np.nan
    return rho, z, zc


def main():
    if not FRAMES:
        print("no frames — run gen_frames.jl first"); return
    ms = [Model(d) for d in FRAMES]; Vs = [meta_V(d) for d in FRAMES]
    Hkm = ms[0].H / 1e3

    # FIXED scales from the fully-loaded final frame: field symlog vmax = the trench dipole peak (so the
    # loaded-edge corner saturates, not the physics); profile x-range spans both columns.
    mf = ms[-1]
    peaks = [np.nanmax(np.abs(rho_profile(mf, xf(mf))[0])) for _, xf, _, _ in COLS]
    vmax = 1.05 * max(peaks)
    rmax = 1.15 * max(peaks)
    norm = SymLogNorm(linthresh=LINTHRESH, vmin=-vmax, vmax=vmax, base=10)
    uzf = mf.array("u", 1); ylo = (-(mf.z[None, :] + WARP * uzf) / 1e3).min() - 3

    for k, (m, V) in enumerate(zip(ms, Vs)):
        ux, uz = m.array("u", 0), m.array("u", 1)
        Xd = (m.x[:, None] + ux) / 1e3                   # horizontal displacement at TRUE scale (as the hero figures); only u_z is exaggerated
        Yd = -(m.z[None, :] + WARP * uz) / 1e3
        C = m.array("dsxz_dx [Pa/m]") / GRAV
        sel = m.x / 1e3 <= WINDOW_KM
        # upsample the SCALAR field (and the deformed coords) before colour-mapping: the blue↔red transition
        # then passes through near-zero (white) cells rather than being colour-blended to grey (gouraud's flaw
        # on the coarse mesh).  order=1 = bounded linear interp, no overshoot.  Smooth AND no grey.
        Xu, Yu, Cu = (zoom(A[sel], (2, 6), order=1) for A in (Xd, Yd, C))
        try:
            pull = trench_pull(m)[0] / 1e12
        except Exception:
            pull = np.nan

        fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.4), gridspec_kw={"width_ratios": [2.15, 1]},
                                       constrained_layout=True)
        # (left) rho_hat field, fixed symlog norm
        axL.pcolormesh(Xu, Yu, Cu, cmap="seismic", norm=norm, shading="gouraud")   # gouraud on the upsampled field: smooth, transitions through white (no grey)
        axL.axhline(0, color="0.5", lw=0.6); axL.set_xlim(0, WINDOW_KM); axL.set_ylim(ylo, 3.0)
        axL.set_xlabel("distance from trench  [km]"); axL.set_ylabel("height  [km]  (deflection ×%g)" % WARP)
        axL.set_title(rf"equivalent density  $\hat\rho=g^{{-1}}\tau_{{zx,x}}$  ·  $V={V:.2f}$ TN m$^{{-1}}$", fontsize=11)
        cb = fig.colorbar(ScalarMappable(norm, "seismic"), ax=axL, pad=0.01, fraction=0.045)
        cb.set_label(r"$\hat\rho$  [kg m$^{-3}$]  (symlog)", fontsize=9)

        # (right) rho_hat(z) at the trench and max-moment columns, with centroids
        for lab, xf, col, lw in COLS:
            rho, z, zc = rho_profile(m, xf(m))
            axR.plot(rho, z, color=col, lw=lw, label=lab)
            if np.isfinite(zc):
                axR.plot(0.0, zc, "o", color=col, ms=9, mec="k", mew=0.9, zorder=6)
        axR.axhline(Hkm / 2, color="0.35", lw=1.0, ls=(0, (5, 3)), alpha=0.8, label="$h/2$")
        axR.axvline(0, color="0.6", lw=0.7); axR.set_xlim(-rmax * 0.55, rmax); axR.set_ylim(Hkm + 2, -2)
        axR.set_xlabel(r"$\hat\rho$  [kg m$^{-3}$]"); axR.set_ylabel("depth below plate surface  [km]")
        axR.set_title("equivalent-density dipole", fontsize=11); axR.legend(fontsize=9, loc="lower right"); axR.grid(alpha=0.25)
        axR.annotate(rf"pull $={pull:.2f}$ TN m$^{{-1}}$", (0.04, 0.05), xycoords="axes fraction",
                     fontsize=10, fontweight="bold")

        p = os.path.join(PNGDIR, f"r_{k:02d}.png"); fig.savefig(p, dpi=110); plt.close(fig)
        print(f"  frame {k:02d}: V={V:.2f}  pull={pull:.2f}")

    pngs = sorted(glob.glob(os.path.join(PNGDIR, "r_*.png")))
    imgs = [Image.open(p).convert("RGB") for p in pngs]
    durs = [140] * len(imgs); durs[-1] = 1400
    imgs[0].save(OUT_GIF, save_all=True, append_images=imgs[1:], duration=durs, loop=0)
    print("wrote", OUT_GIF)
    try:
        import imageio.v2 as iio
        iio.mimsave(OUT_MP4, [iio.imread(p) for p in pngs], fps=8); print("wrote", OUT_MP4)
    except Exception as e:
        print("mp4 skipped:", e)


if __name__ == "__main__":
    main()
