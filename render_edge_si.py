"""render_edge_si.py — SUPPLEMENTARY figure: the near-trench edge effect, why we exclude a thin (<10 km) layer,
and what excluding it does to the compensation arm / trench pull.

(a) pseudo-density field near the trench (Tresca ref) — the perturbation is a thin skin in the top few km at the
    loaded corner, with the moment arm overlaid.
(b) the arm rises through that ~8 km edge window then plateaus; we exclude the window and read the plateau (the
    slow fall to h/2 further inboard is genuine flexure, not edge).
(c) edge-excluded arm ≈ 0.55 h across all four rheologies; excluding the edge shifts the arm only a few %.

    /opt/anaconda3/envs/pyvista-env/bin/python render_edge_si.py
"""
import sys; sys.path.insert(0, ".")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from gpe_analysis import Model, trench_ref_km, deformed_shear_gradient, reference_lines
plt.rcParams.update({"font.size": 11, "axes.titlesize": 11.5, "axes.labelsize": 11})
G, DR, H, THR = 9.81, 2300.0, 60.0, 0.15
S1 = "out/paper/set1_rheology/"
REF = S1 + "tresca_deep_150_60km_V4"
SUITE = [("elastic", S1 + "elastic_deep_60km_V4", "#7f7f7f"),
         ("Tresca", S1 + "tresca_deep_150_60km_V4", "#1a1a1a"),
         ("DD-VM\nasym", S1 + "dd_vm_asym_60km_V4", "#2ca02c"),
         ("DD-VM\nsym", S1 + "dd_vm_sym_60km_V4", "#9467bd")]

def profile(p):
    # arm = DIRECT ΔGPE*(x)/(gΔρw(x)) (self-consistent with the certified pull); skin (ρ̂) only sets L_edge
    m = Model(p); xt = trench_ref_km(m); wT = float(np.nanmax(m.topography()))
    xi = reference_lines(m)["isostatic"]; SzzI = m.deformed_resultants(xi)[0]; top = m.topography()
    XS = np.linspace(max(xt, 0.4), 90, 200); arm, skin = [], []
    for x in XS:
        L, tau = deformed_shear_gradient(m, x); d = (L["z"] - L["z"][0]) / 1e3; skin.append(np.mean(tau[d <= 8.0]) / G)
        Szz = m.deformed_resultants(x)[0]; w = np.interp(x, m.xkm, top)
        arm.append((Szz - SzzI) / (G * DR * w) / 1e3 if w > 50 else np.nan)     # DIRECT arm [km]
    arm, skin = np.array(arm), np.array(skin)
    sfar = np.median(skin[(XS >= 40) & (XS <= 80)])
    dec = np.where((XS > xt) & (np.abs(skin - sfar) < THR * abs(skin[0] - sfar)))[0]
    Le = float(XS[dec[0]]) if dec.size else 10.0
    win = (XS >= Le) & (XS <= 2 * Le)
    return dict(m=m, xt=xt, XS=XS, arm=arm, Le=Le, plat=np.nanmedian(arm[win]),
                lo=np.nanmin(arm[win]), hi=np.nanmax(arm[win]), raw=arm[0], w=wT)

ref = profile(REF)
suite = [(lab, profile(p), c) for lab, p, c in SUITE]

fig = plt.figure(figsize=(11.5, 8.4))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.92], hspace=0.32, wspace=0.26,
                       left=0.08, right=0.95, top=0.93, bottom=0.09)
# (a) field
axa = fig.add_subplot(gs[0, 0]); m = ref["m"]
rho = (m.array("dsxz_dx [Pa/m]") if m.has_field("dsxz_dx [Pa/m]") else m.grad_x_field(m.array("sigma_xz [Pa]"))) / G
xkm, dep = m.xkm, m.depth; sel = xkm <= 40
cl = np.nanpercentile(np.abs(rho), 97)
pc = axa.pcolormesh(xkm[sel], dep, rho[sel].T, cmap="RdBu_r", vmin=-cl, vmax=cl, shading="auto")
a40 = ref["XS"] <= 40; axa.plot(ref["XS"][a40], ref["arm"][a40], "k-", lw=2.4)
axa.add_patch(plt.Rectangle((0, 0), ref["Le"], 10, fill=False, ec="k", lw=1.4, ls=(0, (4, 2))))
axa.annotate("edge skin: thin ($\\lesssim$10 km)\nlayer at the loaded corner", (ref["Le"], 10), (14, 20),
             fontsize=9.5, arrowprops=dict(arrowstyle="->", lw=1.1))
axa.axhline(H / 2, color="0.3", ls="--", lw=1)
axa.invert_yaxis(); axa.set_xlabel("distance from trench [km]"); axa.set_ylabel("depth below top [km]")
axa.set_title("(a)  pseudo-density $\\hat{\\rho}=\\partial_x\\sigma_{xz}/g$  +  moment arm", loc="left")
fig.colorbar(pc, ax=axa, label="$\\hat{\\rho}$ [kg m$^{-3}$]", shrink=.85, pad=.02)
# (b) exclude + plateau
axb = fig.add_subplot(gs[0, 1])
axb.plot(ref["XS"], ref["arm"], "k-", lw=2.6)
axb.axvspan(0, ref["Le"], color="#b04a2f", alpha=.15)
axb.text(ref["Le"] / 2, 36.4, "excluded\nedge", ha="center", fontsize=9.5, color="#8a3320", fontweight="bold")
axb.axhline(ref["plat"], color="#2e7a62", ls="--", lw=1.6)
axb.text(20, ref["plat"] + .22, "reported arm = %.0f km (%.2f h)" % (ref["plat"], ref["plat"] / H), color="#2e7a62", fontsize=9.5)
axb.axhline(H / 2, color="0.55", ls=":", lw=1.4); axb.text(82, H / 2, "h/2", color="0.5", va="center", fontsize=10)
axb.annotate("slow fall to h/2 =\ngenuine flexure, not edge", (60, 32), (40, 34.6), fontsize=9,
             color="0.35", arrowprops=dict(arrowstyle="->", color="0.5"))
axb.set_xlim(0, 84); axb.set_ylim(29, 37); axb.set_xlabel("distance from trench [km]"); axb.set_ylabel("moment arm [km]")
axb.set_title("(b)  exclude the $\\lesssim$8 km edge window, read the plateau", loc="left")
# (c) across rheology
axc = fig.add_subplot(gs[1, :])
x = np.arange(len(suite))
axc.axhspan(0.50, 0.60, color="#2e7a62", alpha=.08)
axc.axhline(0.55, color="#2e7a62", lw=1, ls="--")
for xi, (lab, r, c) in zip(x, suite):
    axc.bar(xi, r["plat"] / H, width=.5, color=c, zorder=3)
    axc.plot(xi, r["raw"] / H, marker="_", ms=26, mew=3, color="0.25", zorder=4)
    axc.text(xi, r["plat"] / H + .006, "%.2f h" % (r["plat"] / H), ha="center", fontsize=10.5, fontweight="bold")
    axc.text(xi, .02, "pull\n%.2f" % (G * DR * r["w"] * r["plat"] / 1e9), ha="center", va="bottom", fontsize=8.5, color="white")
axc.plot([], [], marker="_", ms=16, mew=3, color="0.25", ls="none", label="raw trench arm (edge included)")
axc.plot([], [], color="#2e7a62", ls="--", label="0.55 h (manuscript coefficient)")
axc.set_xticks(x); axc.set_xticklabels([s[0] for s in suite]); axc.set_ylim(0, 0.66)
axc.set_ylabel("edge-excluded arm  /  h"); axc.legend(fontsize=9.5, loc="upper right")
axc.set_title("(c)  edge-excluded arm ≈ 0.55 h across rheology; excluding the edge (bar vs $-$) shifts it only a few %  ·  pull in TN/m", loc="left")
fig.suptitle("Supplementary: the near-trench edge effect is a thin ($\\lesssim$10 km) excludable layer; the compensation arm is ~0.55 h", fontsize=12.5)
fig.savefig("out/paper/figures/edge_effect_SI.png", dpi=200)
print("wrote out/paper/figures/edge_effect_SI.png")
print("\nSuggested caption:")
print("Near-trench edge effect. (a) The pseudo-density field ρ̂=∂ₓσ_xz/g near the loaded trench (Tresca reference,")
print("h=60 km): the perturbation is a thin skin confined to the top ≲10 km at the corner (dashed box); the black")
print("line is the signed moment arm. (b) The arm rises through the ~8 km edge window (identified from the decay of")
print("the top-skin amplitude) and then plateaus; we exclude the window and report the plateau. The slower fall")
print("toward h/2 further inboard is genuine flexural variation, not an edge artefact. (c) The edge-excluded arm is")
print("0.51–0.60 h (≈0.55 h) across all four rheologies; excluding the edge (bar vs dash) changes it by only a few")
print("percent, and the residual spread is rheological. This supports the coefficient 0.55±0.10 used in the main text.")
