"""animation/render_anim_sxz.py — third of the trio (STANDALONE, DISPOSABLE).

Same 24 load-step frames (animation/frames/step_NN), no re-solving.  Shows the VERTICAL SHEAR STRESS
sigma_xz — the parent of the shear-stress gradient — through the loading cycle:
  (left)  sigma_xz as a FIELD on the deformed plate (the bending shear, growing with load);
  (right) sigma_xz(z) DEPTH PROFILES at the TRENCH (single-signed: integrates to the applied shear V) and
          the MAX-BENDING-MOMENT column (S-shaped: integrates to ~0, since V = dM/dx = 0 there).

Scalar upsampled before colour-mapping so blue<->red transitions pass through white (no grey), as in
render_anim_rho.py.  Reuses ../gpe_analysis.py READ-ONLY.  Delete animation/ to remove everything.

  python animation/render_anim_sxz.py
"""
import sys, os, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.ndimage import zoom
from PIL import Image
from gpe_analysis import Model, trench_pull, trench_ref_km, reference_lines

HERE = os.path.dirname(__file__)
FRAMES = sorted(glob.glob(os.path.join(HERE, "frames", "step_*")))
PNGDIR = os.path.join(HERE, "frames", "_png_sxz"); os.makedirs(PNGDIR, exist_ok=True)
OUT_GIF = os.path.join(HERE, "trench_pull_sxz.gif")
OUT_MP4 = os.path.join(HERE, "trench_pull_sxz.mp4")
WINDOW_KM = 260
WARP = 6.0
COLS = [("trench",  trench_ref_km, "k",       2.6),
        ("max $M$", lambda m: reference_lines(m)["moment_max"], "#c0392b", 2.0)]


def meta_V(d):
    try:
        return float(open(os.path.join(d, "meta.txt")).read().split("V_TN=")[1].split()[0])
    except Exception:
        return np.nan


def main():
    if not FRAMES:
        print("no frames — run gen_frames.jl first"); return
    ms = [Model(d) for d in FRAMES]; Vs = [meta_V(d) for d in FRAMES]
    Hkm = ms[0].H / 1e3

    # FIXED linear scales from the fully-loaded final frame
    mf = ms[-1]
    clim = 0.9 * np.nanpercentile(np.abs(mf.array("sigma_xz [Pa]") / 1e6), 98)
    smax = 1.12 * max(np.nanmax(np.abs(mf.deformed_line(xf(mf))["sxz"])) for _, xf, _, _ in COLS) / 1e6
    uzf = mf.array("u", 1); ylo = (-(mf.z[None, :] + WARP * uzf) / 1e3).min() - 3
    norm = Normalize(-clim, clim)

    for k, (m, V) in enumerate(zip(ms, Vs)):
        ux, uz = m.array("u", 0), m.array("u", 1)
        Xd = (m.x[:, None] + WARP * ux) / 1e3
        Yd = -(m.z[None, :] + WARP * uz) / 1e3
        C = m.array("sigma_xz [Pa]") / 1e6
        sel = m.x / 1e3 <= WINDOW_KM
        Xu, Yu, Cu = (zoom(A[sel], (2, 6), order=1) for A in (Xd, Yd, C))    # upsample scalar before colour-map (no grey)
        try:
            pull = trench_pull(m)[0] / 1e12
        except Exception:
            pull = np.nan

        fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.4), gridspec_kw={"width_ratios": [2.15, 1]},
                                       constrained_layout=True)
        axL.pcolormesh(Xu, Yu, Cu, cmap="PuOr_r", norm=norm, shading="gouraud")
        axL.axhline(0, color="0.5", lw=0.6); axL.set_xlim(0, WINDOW_KM); axL.set_ylim(ylo, 3.0)
        axL.set_xlabel("distance from trench  [km]"); axL.set_ylabel("height  [km]  (deflection ×%g)" % WARP)
        axL.set_title(rf"vertical shear stress  $\sigma_{{xz}}$  ·  $V={V:.2f}$ TN m$^{{-1}}$", fontsize=11)
        cb = fig.colorbar(ScalarMappable(norm, "PuOr_r"), ax=axL, pad=0.01, fraction=0.045)
        cb.set_label(r"$\sigma_{xz}$  [MPa]", fontsize=9)

        for lab, xf, col, lw in COLS:
            L = m.deformed_line(xf(m))
            axR.plot(L["sxz"] / 1e6, (L["z"] - L["z"][0]) / 1e3, color=col, lw=lw, label=lab)
        axR.axhline(Hkm / 2, color="0.35", lw=1.0, ls=(0, (5, 3)), alpha=0.8, label="$h/2$")
        axR.axvline(0, color="0.6", lw=0.7); axR.set_xlim(-smax, smax); axR.set_ylim(Hkm + 2, -2)
        axR.set_xlabel(r"$\sigma_{xz}$  [MPa]"); axR.set_ylabel("depth below plate surface  [km]")
        axR.set_title("vertical shear-stress profile", fontsize=11); axR.legend(fontsize=9, loc="lower right"); axR.grid(alpha=0.25)
        axR.annotate(rf"pull $={pull:.2f}$ TN m$^{{-1}}$", (0.04, 0.05), xycoords="axes fraction",
                     fontsize=10, fontweight="bold")

        p = os.path.join(PNGDIR, f"s_{k:02d}.png"); fig.savefig(p, dpi=110); plt.close(fig)
        print(f"  frame {k:02d}: V={V:.2f}  pull={pull:.2f}")

    pngs = sorted(glob.glob(os.path.join(PNGDIR, "s_*.png")))
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
