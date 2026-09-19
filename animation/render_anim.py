"""animation/render_anim.py — STANDALONE, DISPOSABLE trench-pull LOADING animation.

Reads the per-load-step states written by gen_frames.jl (animation/frames/step_NN/gpe_model.vtu) and
builds a 2-panel movie of the loading cycle:
  (top)    the deformed plate (vertical warp), coloured by the normal-stress difference sigma_xx - sigma_zz,
           with the yield front (black contour) growing inward as the elastic core collapses;
  (bottom) the trench deflection w_T and the trench pull dGPE* building up as the end load V ramps.

Dependencies: matplotlib + PIL (both standard) and ../gpe_analysis.py, imported READ-ONLY for the
deformed-Cauchy field extraction.  It modifies nothing in the core; delete animation/ to remove it all.

  python animation/render_anim.py
"""
import sys, os, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))   # ../gpe_analysis.py (read-only)
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from PIL import Image
from gpe_analysis import Model, trench_pull

HERE = os.path.dirname(__file__)
FRAMES = sorted(glob.glob(os.path.join(HERE, "frames", "step_*")))
PNGDIR = os.path.join(HERE, "frames", "_png"); os.makedirs(PNGDIR, exist_ok=True)
OUT_GIF = os.path.join(HERE, "trench_pull_loading.gif")
OUT_MP4 = os.path.join(HERE, "trench_pull_loading.mp4")
WINDOW_KM = 300
WARP = 6.0                                        # vertical deflection exaggeration for the movie
GRAV = 9.81


def load_all():
    rows = []
    for d in FRAMES:
        m = Model(d)
        try:
            V = float(open(os.path.join(d, "meta.txt")).read().split("V_TN=")[1].split()[0])
        except Exception:
            V = np.nan
        w = m.topography().max() / 1e3
        try:
            pull = trench_pull(m)[0] / 1e12
        except Exception:
            pull = np.nan
        rows.append(dict(m=m, V=V, w=w, pull=pull))
    return rows


def main():
    if not FRAMES:
        print("no frames — run gen_frames.jl first"); return
    R = load_all()
    Vs = np.array([r["V"] for r in R]); ws = np.array([r["w"] for r in R]); pulls = np.array([r["pull"] for r in R])

    # fixed colour scale + axes from the final (fully loaded) frame
    mf = R[-1]["m"]
    dif = (mf.array("sigma_xx [Pa]") - mf.array("sigma_zz [Pa]")) / 1e6
    clim = 0.85 * np.nanpercentile(np.abs(dif), 98)
    uzf = mf.array("u", 1); ydef_f = -(mf.z[None, :] + WARP * uzf) / 1e3
    ylo, yhi = ydef_f.min() - 3, 3.0

    for k, r in enumerate(R):
        m = r["m"]
        ux, uz = m.array("u", 0), m.array("u", 1)
        Xd = (m.x[:, None] + ux) / 1e3                   # horizontal displacement at TRUE scale (as the hero figures); only u_z is exaggerated
        Yd = -(m.z[None, :] + WARP * uz) / 1e3
        C = (m.array("sigma_xx [Pa]") - m.array("sigma_zz [Pa]")) / 1e6
        yld = m.array("yielded")
        sel = m.x / 1e3 <= WINDOW_KM

        fig, (axT, axB) = plt.subplots(2, 1, figsize=(9.5, 7.2),
                                       gridspec_kw={"height_ratios": [2.05, 1]}, constrained_layout=True)
        axT.pcolormesh(Xd[sel], Yd[sel], C[sel], cmap="BrBG", vmin=-clim, vmax=clim, shading="gouraud")
        if np.nanmax(yld) > 0.5:
            axT.contour(Xd[sel], Yd[sel], yld[sel], [0.5], colors="k", linewidths=1.1)
        axT.axhline(0, color="0.5", lw=0.6)
        axT.set_xlim(0, WINDOW_KM); axT.set_ylim(ylo, yhi)
        axT.set_ylabel("height  [km]  (deflection ×%g)" % WARP); axT.set_xticklabels([])
        axT.set_title(rf"Trench-pull loading  ·  $V={r['V']:.2f}$ TN m$^{{-1}}$  ·  $\sigma_{{xx}}-\sigma_{{zz}}$ (yield front in black)", fontsize=11)
        cb = fig.colorbar(ScalarMappable(Normalize(-clim, clim), "BrBG"), ax=axT, pad=0.01, fraction=0.05)
        cb.set_label(r"$\sigma_{xx}-\sigma_{zz}$  [MPa]", fontsize=9)

        # bottom: pull (left) and deflection (right) accumulating vs load V
        axB.plot(Vs[:k + 1], pulls[:k + 1], "-o", color="#c0392b", ms=4, lw=1.8, label=r"pull $\Delta\mathrm{GPE}^{*}$")
        axB.plot(Vs[k], pulls[k], "o", color="#c0392b", ms=9, mec="k", mew=0.8, zorder=5)
        axB.set_xlim(0, Vs.max() * 1.05); axB.set_ylim(0, np.nanmax(pulls) * 1.12)
        axB.set_xlabel(r"applied load  $V$  [TN m$^{-1}$]"); axB.set_ylabel(r"$\Delta\mathrm{GPE}^{*}$  [TN m$^{-1}$]", color="#c0392b")
        axB.tick_params(axis="y", labelcolor="#c0392b"); axB.grid(alpha=0.25)
        axR = axB.twinx()
        axR.plot(Vs[:k + 1], ws[:k + 1], "-s", color="#1f3b73", ms=3.5, lw=1.5)
        axR.plot(Vs[k], ws[k], "s", color="#1f3b73", ms=8, mec="k", mew=0.7, zorder=5)
        axR.set_ylim(0, np.nanmax(ws) * 1.12); axR.set_ylabel(r"trench deflection  $w_T$  [km]", color="#1f3b73")
        axR.tick_params(axis="y", labelcolor="#1f3b73")
        axB.annotate(rf"$w_T={r['w']:.2f}$ km   pull $={r['pull']:.2f}$ TN m$^{{-1}}$",
                     (0.03, 0.9), xycoords="axes fraction", fontsize=10, fontweight="bold")

        p = os.path.join(PNGDIR, f"f_{k:02d}.png"); fig.savefig(p, dpi=110); plt.close(fig)
        print(f"  frame {k:02d}: V={r['V']:.2f}  w={r['w']:.2f}km  pull={r['pull']:.2f}")

    pngs = sorted(glob.glob(os.path.join(PNGDIR, "f_*.png")))
    imgs = [Image.open(p).convert("RGB") for p in pngs]
    durs = [140] * len(imgs); durs[-1] = 1400                          # hold on the fully-loaded frame
    imgs[0].save(OUT_GIF, save_all=True, append_images=imgs[1:], duration=durs, loop=0)
    print("wrote", OUT_GIF)
    try:                                                              # mp4 if imageio+ffmpeg present
        import imageio.v2 as iio
        iio.mimsave(OUT_MP4, [iio.imread(p) for p in pngs], fps=8)
        print("wrote", OUT_MP4)
    except Exception as e:
        print("mp4 skipped (no imageio-ffmpeg):", e)


if __name__ == "__main__":
    main()
