"""render_gpe_correlation.py — is the trench pull a purely TOPOGRAPHIC effect?

Correlates ΔGPE* (trench pull) against the applied load V and the trench deflection w, across BOTH the
load sweep (Suite 2, N_D=0, V varied) and the in-plane sweep (Suite 3, V=3.5, N_D varied).  The key test:
if ΔGPE* is purely topographic it collapses onto ONE ΔGPE*-vs-w curve regardless of whether w was produced
by changing V or N_D.

  (A) ΔGPE* vs V        — Suite 2 sweeps V; Suite 3 sits at V=3.5 (so N_D spreads the pull at fixed load).
  (B) ΔGPE* vs w        — both suites; do they collapse onto one line?  (linearity in topography)
  (C) V vs w            — the load→deflection map (yielding bends it away from linear).

    /opt/anaconda3/envs/pyvista-env/bin/python render_gpe_correlation.py
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from gpe_analysis import Model, trench_pull

P = "data"
# (V, N_D, path)
# Baseline is now V=4 (σ_Y=150, M_p=1.35).  V-sweep spans V1..V4.5 in suite2_load; the V4 point IS the baseline.
SUITE2 = [(1.0, 0, f"{P}/suite2_load/tresca_deep_150_60km_V1"), (1.5, 0, f"{P}/suite2_load/tresca_deep_150_60km_V1p5"),
          (2.0, 0, f"{P}/suite2_load/tresca_deep_150_60km_V2"), (2.5, 0, f"{P}/suite2_load/tresca_deep_150_60km_V2p5"),
          (3.0, 0, f"{P}/suite2_load/tresca_deep_150_60km_V3"), (3.5, 0, f"{P}/suite2_load/tresca_deep_150_60km_V3p5"),
          (4.0, 0, f"{P}/suite1_strength/tresca_deep_150_60km_V4"), (4.5, 0, f"{P}/suite2_load/tresca_deep_150_60km_V4p5")]
SUITE3 = [(4.0, -3, f"{P}/suite3_background/tresca_deep_150_60km_V4_mem-3"), (4.0, -2, f"{P}/suite3_background/tresca_deep_150_60km_V4_mem-2"),
          (4.0, -1, f"{P}/suite3_background/tresca_deep_150_60km_V4_mem-1"), (4.0, 1, f"{P}/suite3_background/tresca_deep_150_60km_V4_mem1"),
          (4.0, 2, f"{P}/suite3_background/tresca_deep_150_60km_V4_mem2"), (4.0, 3, f"{P}/suite3_background/tresca_deep_150_60km_V4_mem3")]
# Suite 1 — rheology (V=4, N_D=0, all M_p=1.35): the Tresca baseline is the shared anchor (the V=4 point in
# Suite 2 / centre of Suite 3), so here we plot only the OTHER rheologies varying around it.
SUITE1 = [("elastic", f"{P}/suite1_strength/elastic_deep_60km_V4"), ("DD-VM asym", f"{P}/suite1_strength/dd_vm_asym_60km_V4"),
          ("DD-VM sym", f"{P}/suite1_strength/dd_vm_sym_60km_V4")]
OUT = f"{P}/reference_figures/gpe_correlation.png"


def gather(rows):
    V, ND, w, G = [], [], [], []
    for v, nd, p in rows:
        m = Model(p); V.append(v); ND.append(nd)
        w.append(m.topography().max() / 1e3); G.append(trench_pull(m)[0] / 1e12)
    return map(np.array, (V, ND, w, G))


def gather1(rows):                                   # Suite 1: label + (w, ΔGPE*), all at V=3.5
    lab, w, G = [], [], []
    for name, p in rows:
        m = Model(p); lab.append(name)
        w.append(m.topography().max() / 1e3); G.append(trench_pull(m)[0] / 1e12)
    return lab, np.array(w), np.array(G)


def r2(x, y):
    a, b = np.polyfit(x, y, 1); yh = a * x + b
    return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2), a, b


def smart_labels(fig, a, obst, xs, ys, texts, color, fs=8.5):
    """Place each label in clear space (never over a data point or another label) with a leader to
    its own point.  obst = all data points to avoid (Nx2 in data coords).  Requires a frozen layout."""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float); obst = np.asarray(obst, float)
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    P = a.transData.transform(np.column_stack([xs, ys]))          # label anchor points (pixels)
    O = a.transData.transform(obst)                               # obstacle points (pixels)
    cand = [(r, ang) for r in (28, 42, 58, 76, 96) for ang in (-90, -55, -125, -30, -150, -18, -162, 22, 158, 48, 132)]
    placed = []
    for idx in np.argsort(-ys):                                   # highest-y first
        px, py = P[idx]
        t = a.text(0.5, 0.5, texts[idx], fontsize=fs, transform=a.transAxes)
        bb = t.get_window_extent(rend); tw, th = bb.width, bb.height; t.remove()
        best, bs = None, -1e18
        for r, ang in cand:
            ox, oy = r * np.cos(np.radians(ang)), r * np.sin(np.radians(ang))
            cx, cy = px + ox, py + oy
            box = (cx - tw / 2 - 7, cy - th / 2 - 5, cx + tw / 2 + 7, cy + th / 2 + 5)
            if np.hypot(O[:, 0] - cx, O[:, 1] - cy).min() < tw / 2 + 12:          # clear of every data point
                continue
            if any(not (box[2] < q[0] or box[0] > q[2] or box[3] < q[1] or box[1] > q[3]) for q in placed):
                continue                                                          # clear of already-placed labels
            score = -r + (10 if oy < 0 else 0)                                    # nearest clear spot, prefer below
            if score > bs:
                bs, best = score, (ox, oy, box)
        if best is None:
            ox, oy = -48, -28; cx, cy = px + ox, py + oy
            best = (ox, oy, (cx - tw / 2, cy - th / 2, cx + tw / 2, cy + th / 2))
        ox, oy, box = best; placed.append(box)
        a.annotate(texts[idx], (xs[idx], ys[idx]), xytext=(ox, oy), textcoords="offset points",
                   fontsize=fs, color=color, ha="center", va="center",
                   path_effects=[pe.withStroke(linewidth=2.6, foreground="white")],
                   arrowprops=dict(arrowstyle="-", color=color, lw=0.7, alpha=0.85, shrinkA=1, shrinkB=6))


def main():
    V2, _, w2, G2 = gather(SUITE2)          # Suite 2: load sweep (N_D=0)
    V3, N3, w3, G3 = gather(SUITE3)          # Suite 3: in-plane N_D sweep (V=4)
    lab1, w1, G1 = gather1(SUITE1)           # Suite 1: rheology (V=4, N_D=0)
    C1, C2, C3 = "#2ca02c", "#1f77b4", "#d62728"     # suite 1 / 2 / 3 colours
    iref = int(np.argmin(np.abs(V2 - 4.0))); wref, Gref = w2[iref], G2[iref]      # reference = the V=4 point (member of all suites)
    fig, ax = plt.subplots(1, 2, figsize=(11.8, 5.1), constrained_layout=True)

    # (labels placed at the end, after the layout is frozen, via module-level smart_labels)

    def suites(a, x1, x2, x3, xref, s1=True, z1=6):                              # s1: show rheology suite; z1: its z-order (put low in panel A)
        if s1:
            a.scatter(x1, G1, c=C1, marker="^", s=62, ec="k", lw=0.4, zorder=z1, label="Suite 1 (strength model)")
        a.plot(x2, G2, "-", color=C2, lw=1.6, zorder=3)                          # continuous V-sweep guide line
        # square = symmetric-Tresca (reference) rheology; colour = what is varied.  All of Suite 2, Suite 3
        # and the reference share this one rheology, so all are squares; only Suite 1 (other rheologies) differs.
        a.plot(x2, G2, "s", color=C2, ms=6, mfc="white", mew=1.4, zorder=4, label="Suite 2 (load $V$)")
        a.scatter(x3, G3, c=C3, marker="s", s=58, ec="k", lw=0.4, zorder=6, label="Suite 3 (in-plane $N_D$)")
        a.plot([xref], [Gref], "s", color=C2, ms=6, zorder=10,
               label="reference (Tresca, $V$=4, $N_D$=0)")
        a.grid(alpha=0.25); a.set_ylabel(r"$\Delta$GPE$^{*}$  [TN m$^{-1}$]")

    # (A) ΔGPE* vs V — the N_D models are labelled here
    suites(ax[0], np.full_like(G1, 4.0), V2, np.full_like(G3, 4.0), 4.0, s1=True, z1=1)
    ax[0].set_xlabel(r"applied load  $V$  [TN m$^{-1}$]")
    ax[0].set_title("(A)  pull vs load", fontsize=11); ax[0].legend(fontsize=8, loc="upper left")

    # (B) ΔGPE* vs w — the rheology suite is labelled here; vs the uniform-plate line ΔGPE* = ΔP_T·(H/2) = Δρg·w·(H/2)
    suites(ax[1], w1, w2, w3, wref)
    DRHOG = (3300.0 - 1000.0) * 9.81
    slope = DRHOG * (60e3 / 2) / 1e9                   # TN/m per km  (Δρg·H/2 = 0.677)
    xw = np.array([0.0, max(w2.max(), w3.max()) * 1.03])
    ax[1].plot(xw, slope * xw, "k--", lw=1.4, label=r"uniform plate:  $\Delta P_T\cdot\frac{h}{2}$")
    ax[1].set_xlabel(r"trench deflection  $w$  [km]")
    ax[1].set_title("(B)  pull vs topography", fontsize=11); ax[1].legend(fontsize=8, loc="upper left")

    fig.suptitle(r"Trench pull  $\Delta$GPE$^{*}$", fontsize=12)

    fig.canvas.draw(); fig.set_layout_engine("none")                             # freeze layout so label placement is exact
    obstA = np.vstack([np.column_stack([V2, G2]), np.column_stack([np.full_like(G3, 4.0), G3]),
                       np.column_stack([np.full_like(G1, 4.0), G1]), [[4.0, Gref]]])
    smart_labels(fig, ax[0], obstA, np.full_like(G3, 4.0), G3, [f"$N_D$={int(round(n)):+d}" for n in N3], C3)
    obstB = np.vstack([np.column_stack([w2, G2]), np.column_stack([w3, G3]), np.column_stack([w1, G1]), [[wref, Gref]]])
    smart_labels(fig, ax[1], obstB, w1, G1, lab1, C1)

    fig.savefig(OUT, dpi=140); print("wrote", OUT.split("/")[-1])
    print(f"  reference (Tresca V=4): w={wref:.2f} km, ΔGPE*={Gref:.2f} TN/m ;  uniform-plate slope Δρg·H/2 = {slope:.3f}")
    print("  Suite 1: " + ", ".join(f"{t} (w={x:.2f}, {y:.2f})" for t, x, y in zip(lab1, w1, G1)))


if __name__ == "__main__":
    if len(sys.argv) > 1:        # optional: render to a preview path instead of the committed figure
        OUT = sys.argv[1]
    main()
