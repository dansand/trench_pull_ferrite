"""render_hero.py — modular 4-panel 'hero' figure for the trench region.

PyVista renders each field map on the warped (deflected) plate; matplotlib composites them with
LaTeX colorbars, dashed reference lines, and a line panel for the resultants:
  (a) topography w with the Winkler support (dV/dx)/Δρg
  map panels, selected by HERO_PANELS: normal-stress difference σxx − σzz (+ yield-front contour + principal-stress
  crosses), vertical shear stress τzx (SI layout only), equivalent density ρ̂ = g⁻¹ τzx,x (+ signed centroid track)
  last panel: the resultants V, M, ΔGPE*, ΔN_D
  (d) resultants V, M, F_D vs x (normalised)
Vertical lines mark the flexure reference locations; an extra DASHED line marks DASH_BETWEEN's midpoint.

    python scripts/render_hero.py
"""
import sys, os; sys.path.insert(0, "analysis")
import numpy as np
import pyvista as pv
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.transforms import blended_transform_factory
import matplotlib.patheffects as pe
from gpe_analysis import Model, reference_lines, trench_ref_km, trench_pull, deformed_shear_gradient, write_table, G, DRHOG, CENTROID_RATIO_MIN

# equivalent density gradient g⁻¹ τzx,x: HERO_BRANCH="B" (default) = FE-exported dsxz_dx; "A" = finite-difference grad_x_field
BRANCH = os.environ.get("HERO_BRANCH", "B").upper()
def pgrad(mod):
    if BRANCH == "B":
        if not mod.has_field("dsxz_dx [Pa/m]"):                 # never fall back silently (audit F7)
            raise SystemExit(f"{mod.dir}: no FE-exported dsxz_dx field. The hero figure uses Branch B (the exporter's "
                             "gradient of the Cauchy shear); re-export the model, or set HERO_BRANCH=A explicitly for a "
                             "finite-difference comparison render.")
        return mod.array("dsxz_dx [Pa/m]")
    print("WARNING: HERO_BRANCH=A — finite-difference gradient of the Cauchy shear; comparison only, the paper uses Branch B")
    return mod.grad_x_field(mod.cauchy_fields()[2])

# ===================== CONFIG (iterate here) =====================
MODEL_DIR = "data/suite1_strength/tresca_deep_150_60km_V4"
WINDOW_KM = float(os.environ.get("HERO_WINDOW_KM", 300))    # env override for thinner plates (shorter flexure)
WARP_FACTOR = 5.0                     # VERTICAL displacement exaggeration; horizontal displacement is drawn at true scale (2026-09-18)
VEXAG     = 1.2                       # vertical scale of the warped geometry (1.5 until 2026-09-18; 1.2 gives a wider page aspect)
GRAV      = G
REF_KEYS  = ["trench", "moment_max", "shear_max", "outer_rise"]
REF_LABELS = {"trench": "trench ($x_T$)", "moment_max": "max $M$", "shear_max": "isostatic ($x_I$)", "outer_rise": "forebulge"}
DASH_BETWEEN = ("moment_max", "shear_max")
DASH_LABEL = "mid"
SHOW_MID = False     # DISABLED 2026-08-05 (user): the 'mid' dashed line is the midpoint between max M and the
                     # isostatic column --- redundant/confusing (max M already sits near the geometric
                     # trench->isostatic midpoint), and removed from the profiles figure too.  Set True to revert.
CROSS_STRIDE = (16, 4)
CROSS_LEN_KM = 8.5
CROSS_MIN_FRAC = 0.12
CROSS_COLORS = ("black", "black")      # black on the red/blue field (sign is shown by the field itself)
CROSS_LW = 5.2                         # principal-stress crosses — bold for publication size
CLIM_PCT = 98
SATURATE  = 0.62
RENDER_W  = 1500                      # px width of each PyVista panel render
MAP_W_IN  = 6.6                       # matplotlib map-axes width [in]; height follows the plate aspect
DIFF_DIR  = None                      # if set, fields become (MODEL_DIR − DIFF_DIR); crosses/centroid dropped
OUT = "figures/hero_tresca_deep60.png"
# Map panels: (key, LaTeX colorbar label, cmap).  HERO_PANELS selects which are drawn, by key letter:
#   d = σxx−σzz with the crosses, s = τzx, r = ρ̂.  Fig. 3 (the default invocation) draws "dr" — the shear panel is
#   omitted because everything on it appears elsewhere (its yield contours on the σxx−σzz panel, its gradient on the
#   ρ̂ panel, its profiles in Fig. 6) and the four-panel figure fits a page with its caption (2026-09-18).  The SI
#   hero figures pass HERO_PANELS=dsr for the full five-panel layout.  Panel letters re-run a, b, c, … in order.
ALL_PANELS = {"d": ("diff", r"$\sigma_{xx}-\sigma_{zz}$  [MPa]", "BrBG"),      # σxx−σzz: brown (compression) – teal (tension)
              "s": ("sxz",  r"$\tau_{zx}$  [MPa]", "PRGn"),                    # vertical shear: purple – green
              "r": ("rho",  r"$\hat{\rho}=g^{-1}\tau_{zx,x}$  [kg m$^{-3}$]", "seismic")}  # stress gradient: blue – red
# One colour map per quantity, everywhere (Fig. 3, the SI heroes, the animations): σxx−σzz brown–teal, τzx purple–green,
# ρ̂ blue–red (2026-09-18: ρ̂ is mostly positive and its faint negative field reads best in pale blue; σxx−σzz and τzx use
# both halves of their maps).
PANEL_KEYS = os.environ.get("HERO_PANELS", "dr")
PANELS = [ALL_PANELS[k] for k in PANEL_KEYS]
# ================================================================


def load(model_dir, diff_dir=None):
    """Load a model; if diff_dir is given, the map fields become (model − diff_model), node-wise on
    the (identical) grid, with geometry/crosses/yield taken from the primary model."""
    m = Model(model_dir)
    sxx, szz, sxz = m.cauchy_fields()                        # Cauchy maps (the VTU's sigma_* are 2nd-PK)
    f = dict(sxx=sxx, szz=szz, sxz=sxz, ux=m.array("u", 0), uy=m.array("u", 1))
    diffN, sxzN, rhoN = (sxx - szz), sxz, pgrad(m)
    m2 = None
    if diff_dir:
        m2 = Model(diff_dir)
        sxx2, szz2, sxz2 = m2.cauchy_fields()
        diffN = diffN - (sxx2 - szz2)
        sxzN = sxzN - sxz2
        rhoN = rhoN - pgrad(m2)
    mesh = m.m
    mesh.point_data["diff"] = (diffN / 1e6)[m._ix, m._iz]
    mesh.point_data["sxz"] = (sxzN / 1e6)[m._ix, m._iz]
    mesh.point_data["rho"] = (rhoN / GRAV)[m._ix, m._iz]
    mesh.point_data["yielded"] = m.array("yielded")[m._ix, m._iz]
    # RENDER IN HEIGHT = −depth: the native mesh is z-DOWN (y=depth), but the maps display surface-up.
    # Negate the vertical coordinate AND the vertical warp component together so the plate shows depth
    # increasing downward and deflects DOWN at the trench (consistent with panel (a)).  Both signs flip on purpose:
    # the VTU is z-down (depth) while VTK draws y up; flipping the coordinate alone would put the deflection the
    # wrong way, flipping the warp alone would draw the plate upside down (audit F5).
    mesh.points[:, 1] = -mesh.points[:, 1]
    u = mesh.point_data["u"]
    # Display displacement: horizontal at TRUE scale, vertical × WARP_FACTOR.  Exaggerating u_x as well drew the
    # rotated trench edge 5× too wide (an 8 km sloping face for a 1.7 km tilt) with the "trench" line through its
    # middle; at true scale the edge tilts by its real amount and the trench column x_T is a vertical line tangent
    # to its lower corner, as in the analysis.
    mesh.point_data["u3"] = np.column_stack([u[:, 0], -WARP_FACTOR * u[:, 1], np.zeros(len(u))])
    mesh.point_data["xcoord"] = mesh.points[:, 0]
    x0 = 0.0                                      # reference edge at X=0; window runs X = 0 → WINDOW_KM
    crop = mesh.threshold((x0, WINDOW_KM * 1e3), scalars="xcoord")
    warped = crop.warp_by_vector("u3", factor=1.0).scale([1, VEXAG, 1], inplace=False)
    lines = reference_lines(m)
    lines["trench"] = trench_ref_km(m)            # the trench COLUMN of the analysis (leftmost complete deformed column), not X = 0
    return m, m2, f, warped, lines, x0


def deviatoric_crosses(m, f, x0):
    """(tension, compression) PolyData of headless principal-stress crosses.  Orientation = the in-plane principal
    axes of the Cauchy stress; both arms have the same physical length, ∝ the maximum in-plane shear (σ1−σ2)/2,
    scaled so the largest cross in the figure has arms of CROSS_LEN_KM.
    The crosses are drawn on the DISPLAYED plate, whose geometry is the reference mesh displaced by (u_x, WARP_FACTOR·u_z),
    flipped to height-up and stretched by VEXAG.  So each arm is a physical vector pushed through the SAME map the
    mesh goes through — the display Jacobian D = S·(I + W·∇u)·F⁻¹ with W = diag(1, WARP_FACTOR), S = diag(1, −VEXAG) and
    F = I + ∇u — for its ORIENTATION only, so a surface-parallel principal axis stays parallel to the drawn surface;
    the arm LENGTH is the physical one, so the two arms of a cross are always equal on the page (they both show
    (σ1−σ2)/2).  The two arms are perpendicular in the physical plate, not exactly on the (exaggerated) page."""
    si, sj = CROSS_STRIDE
    ii = [i for i in range(0, m.Nx, si) if m.x[i] >= x0]
    jj = list(range(0, m.Nz, sj))
    ux, uz = f["ux"], f["uy"]
    dux_dx, dux_dz = np.gradient(ux, m.x, axis=0), np.gradient(ux, m.z, axis=1)
    duz_dx, duz_dz = np.gradient(uz, m.x, axis=0), np.gradient(uz, m.z, axis=1)
    S = np.diag([1.0, -VEXAG]); I2 = np.eye(2); Wm = np.diag([1.0, WARP_FACTOR])
    recs, dev = [], []
    for i in ii:
        for j in jj:
            T = np.array([[f["sxx"][i, j], f["sxz"][i, j]], [f["sxz"][i, j], f["szz"][i, j]]])
            wv, V = np.linalg.eigh(T)                                  # physical (x, z-down) principal axes of the Cauchy stress
            Gu = np.array([[dux_dx[i, j], dux_dz[i, j]], [duz_dx[i, j], duz_dz[i, j]]])
            Dmap = S @ (I2 + Wm @ Gu) @ np.linalg.inv(I2 + Gu)       # physical direction → displayed direction
            xw = m.x[i] + ux[i, j]                                     # horizontal displacement at true scale
            yw = -(m.z[j] + WARP_FACTOR * uz[i, j]) * VEXAG            # height = −depth (see load)
            recs.append((xw, yw, V, wv - wv.mean(), Dmap)); dev.append(np.max(np.abs(wv - wv.mean())))
    maxdev = max(dev); scale = (CROSS_LEN_KM * 1e3) / (maxdev + 1e-30)
    seg = {"t": ([], []), "c": ([], [])}
    for xw, yw, V, dv, Dmap in recs:
        if np.max(np.abs(dv)) < CROSS_MIN_FRAC * maxdev:
            continue
        for k in range(2):
            d = Dmap @ V[:, k]; d /= (np.hypot(*d) + 1e-30)            # ORIENTATION follows the drawn plate (display map) …
            a = scale * abs(dv[k]) * d                                 # … LENGTH stays the physical one: both arms equal on the page
            pts, lines = seg["t" if dv[k] > 0 else "c"]
            n = len(pts)
            pts += [[xw - a[0], yw - a[1], 0.0], [xw + a[0], yw + a[1], 0.0]]
            lines += [2, n, n + 1]

    def mk(pl):
        pts, lines = pl
        if not pts:
            return None
        poly = pv.PolyData(); poly.points = np.array(pts); poly.lines = np.array(lines)
        return poly
    return mk(seg["t"]), mk(seg["c"])


def render_field(warped, key, cmap, clim, contour, tcross, ccross, bbox):
    """Render one field map on the warped plate; return an RGB image array framed exactly to bbox."""
    x_ext, y_ext = bbox[1] - bbox[0], bbox[3] - bbox[2]
    asp = x_ext / y_ext
    pl = pv.Plotter(off_screen=True, window_size=(RENDER_W, int(RENDER_W / asp)), border=False)
    pl.background_color = "white"
    pl.add_mesh(warped, scalars=key, cmap=cmap, clim=clim, show_scalar_bar=False)
    if contour is not None and contour.n_points:
        pl.add_mesh(contour, color="black", line_width=5)
    if key == "diff":                      # principal-stress crosses only on (a)
        if tcross is not None:
            pl.add_mesh(tcross, color=CROSS_COLORS[0], line_width=CROSS_LW)
        if ccross is not None:
            pl.add_mesh(ccross, color=CROSS_COLORS[1], line_width=CROSS_LW)
    pl.enable_parallel_projection(); pl.view_xy()
    pl.camera.focal_point = (0.5 * (bbox[0] + bbox[1]), 0.5 * (bbox[2] + bbox[3]), 0.0)
    pl.camera.parallel_scale = y_ext / 2
    img = pl.screenshot(return_img=True); pl.close()
    return img


def resultants(m):
    """Panel-(e) resultants on the DEFORMED column (Cauchy, Model.deformed_line), trench-referenced — the
    one definition.  V=∫σxz, ΔN_D=∫(σxx−σzz)−(trench), M=∫σxx(z−z_mid) (PURE bending moment, so dM/dx=V holds
    exactly — a moment of σxx−σzz couples in σzz and breaks it), ΔGPE*=−(∫σzz−∫σzz_trench)."""
    x_tr = trench_ref_km(m)
    xs = np.linspace(x_tr, WINDOW_KM, 90)
    Szz0 = m.deformed_resultants(x_tr)[0]
    V, ND, M, dG = [], [], [], []
    for xk in xs:
        L = m.deformed_line(xk); z = L["z"]; zmid = 0.5 * (z[0] + z[-1])
        V.append(np.trapz(L["sxz"], z))
        ND.append(np.trapz(L["n_d"], z))
        M.append(np.trapz(L["sxx"] * (z - zmid), z))          # PURE bending moment ∫σxx(z−z_mid) ⇒ dM/dx = V exactly
        dG.append(-(np.trapz(L["szz"], z) - Szz0))
    ND = np.array(ND); ND = ND - ND[0]                         # ΔN_D, trench-referenced
    return xs, np.array(V), ND, np.array(M), np.array(dG)


def build(model_dir=None, out=None, diff_dir=None, save=True):
    """Draw the hero figure for `model_dir` (default: the module's MODEL_DIR). Returns the matplotlib
    figure; writes it to `out` (default OUT) when save=True. From a notebook: build(model_dir, save=False)."""
    model_dir = MODEL_DIR if model_dir is None else model_dir
    diff_dir = DIFF_DIR if diff_dir is None else diff_dir
    out_path = OUT if out is None else out
    m, m2, f, warped, lines, x0 = load(model_dir, diff_dir)
    is_diff = m2 is not None
    tcross, ccross = (None, None) if is_diff else deviatoric_crosses(m, f, x0)   # crosses meaningless for a difference
    contour = warped.contour([0.5], scalars="yielded") if np.nanmax(warped.point_data["yielded"]) > 0.5 else None
    bbox = warped.bounds
    print("reference lines [km]:", {k: round(lines[k], 0) for k in REF_KEYS})

    # SIGNED ρ̂ centroid = the dipole-moment depth of τ_zx,x, via the deformed-line extractor (the same
    # estimator as thickness_compare): per column τ = g⁻¹ τzx,x on the DEFORMED line, signed centroid depth
    # z_c = ∫zτ dz / ∫τ dz.  Only COMPLETE columns count (the track therefore starts at the trench column, ≈1 km in).  Mask ONLY where the dipole is near-BALANCED — |∫τ| < RATIO_MIN·∫|τ| — the genuine
    # ill-conditioning of a signed centroid (the isostatic column, and a fully-plastic hinge where the ± lobes
    # nearly cancel).  This replaces the old reference-grid column-sum ρ̂.sum(depth), which is UNFAITHFUL at a
    # strongly-rotated yielded hinge: on the h=30 plate it flipped sign where the true net resultant is finite
    # and positive, planting spurious zero-crossings that blanked the whole hinge.
    # The hero track is a CONTINUOUS centroid line over the window, masked column by column on x-smoothed column
    # integrals; it keeps its own guard (RATIO_MIN below) so the line runs as far as the dipole is meaningfully signed.
    # The single-column MARKERS in the other figures use the shared CENTROID_RATIO_MIN (0.45) from gpe_analysis on
    # raw columns.  Changing these values moves the ends of the track in Figs 3 and S6–S8 — re-approve those figures.
    RATIO_MIN = 0.25                                       # was 0.35 on raw columns; 0.25 on the SMOOTHED integrals (below)
    SMOOTH_KM = 6.0                                        # Gaussian σ along x applied to the column integrals Q, M1, A before
                                                           # dividing (3 columns): damps the mesh-scale ripple of the deformed-line
                                                           # extraction so the track runs closer to the isostatic column (gap 30 → 20 km
                                                           # on the reference plate) and the 30 km hinge masks as one clean gap.  The
                                                           # pole itself (Q = 0 at x_I) is intrinsic and stays masked.  (2026-09-18)
    wsel_c = m.xkm <= WINDOW_KM
    zc_depth = np.full(m.Nx, np.nan)                       # centroid depth below the deformed top [km]
    yc_warp = np.full(m.Nx, np.nan)                        # centroid warped height [m] (for the map overlay)
    Qc = np.full(m.Nx, np.nan); M1c = np.full(m.Nx, np.nan); Ac = np.full(m.Nx, np.nan)
    zline = {}
    for i in np.where(wsel_c)[0]:
        try:
            L, tau = deformed_shear_gradient(m, m.xkm[i])
        except Exception:
            continue
        z = L["z"]
        if z.size < m.Nz:                                  # the line must cut the WHOLE plate: a vertical line at
            continue                                       # reference x = 0 exits through the rotated edge face part way
                                                           # down (23 of 49 layers) and its "centroid" is that of the
                                                           # truncated column — the false shallowing seen before 2026-09-18
        Qc[i] = np.trapz(tau, z); Ac[i] = np.trapz(np.abs(tau), z); M1c[i] = np.trapz(z * tau, z); zline[i] = L
    ok = np.isfinite(Qc)
    if ok.sum() > 3:
        from scipy.ndimage import gaussian_filter1d
        sig = SMOOTH_KM / (m.xkm[1] - m.xkm[0])
        Qc[ok] = gaussian_filter1d(Qc[ok], sig); M1c[ok] = gaussian_filter1d(M1c[ok], sig); Ac[ok] = gaussian_filter1d(Ac[ok], sig)
    for i in np.where(ok)[0]:
        Q, A = Qc[i], Ac[i]
        if A <= 0 or abs(Q) < RATIO_MIN * A:              # near-balanced ⇒ signed centroid ill-conditioned ⇒ skip
            continue
        zc = M1c[i] / Q                                    # absolute deformed depth of the centroid
        L = zline[i]; z = L["z"]
        if not (z[0] < zc < z[-1]):
            continue
        zref_c = np.interp(zc, z, L["z_ref"]); uy_c = zc - zref_c   # map deformed depth → (reference depth, u_y)
        zc_depth[i] = (zc - z[0]) / 1e3
        yc_warp[i] = -(zref_c + WARP_FACTOR * uy_c) * VEXAG          # warped height, consistent with the mesh warp
    print(f"rho-hat centroid depth over window: mean {np.nanmean(zc_depth[wsel_c]):.1f} km "
          f"(range {np.nanmin(zc_depth[wsel_c]):.1f}-{np.nanmax(zc_depth[wsel_c]):.1f})")

    # --- render the three field maps with PyVista ---
    imgs, clims = [], []
    for key, _, cmap in PANELS:
        v = warped.point_data[key]
        cl = SATURATE * np.nanpercentile(np.abs(v), CLIM_PCT)
        imgs.append(render_field(warped, key, cmap, (-cl, cl), contour, tcross, ccross, bbox))
        clims.append(cl)

    # --- composite in matplotlib ---
    x0km, Lkm = bbox[0] / 1e3, bbox[1] / 1e3
    y0km, y1km = bbox[2] / 1e3, bbox[3] / 1e3
    aspect = (bbox[1] - bbox[0]) / (bbox[3] - bbox[2])
    map_h = MAP_W_IN / aspect
    line_h = map_h                                    # all five panels the same height (2026-09-18)
    top_h = map_h
    nmap = len(PANELS)
    fig = plt.figure(figsize=(MAP_W_IN + 1.5, top_h + nmap * map_h + line_h + 1.2), constrained_layout=True)
    gs = fig.add_gridspec(nmap + 2, 2, width_ratios=[MAP_W_IN, 0.22],
                          height_ratios=[top_h] + [map_h] * nmap + [line_h])

    x_dash = 0.5 * (lines[DASH_BETWEEN[0]] + lines[DASH_BETWEEN[1]])
    refs = [(lines[k], REF_LABELS[k]) for k in REF_KEYS if x0km <= lines[k] <= Lkm]
    wsel = (m.xkm >= x0km) & (m.xkm <= Lkm)
    lam = 2 * np.pi * ((4 * (70e9 / (1 - 0.25**2)) * (m.H**3 / 12.0) / DRHOG) ** 0.25)  # flexural wavelength 2πα [m]

    def verticals(ax, top_labels=False):
        for xk, lab in refs:
            ax.axvline(xk, color="0.25", lw=1.4)
        if SHOW_MID:
            ax.axvline(x_dash, color="0.25", lw=1.9, ls=(0, (7, 4)))
        if top_labels:
            tr = blended_transform_factory(ax.transData, ax.transAxes)
            for xk, lab in refs:
                ax.text(xk, 1.02, lab, transform=tr, ha="center", va="bottom", fontsize=12)
            if SHOW_MID:
                ax.text(x_dash, 1.02, DASH_LABEL, transform=tr, ha="center", va="bottom", fontsize=12, color="0.35")

    # panel (a): topography and the foundation support (dV/dx)/(Δρg) — these coincide
    axt = fig.add_subplot(gs[0, 0])
    if is_diff:
        wtopo, supp = m.topography() - m2.topography(), m.grad_x(m.V() - m2.V()) / DRHOG
        lw_, ls_, ylt = r"$\Delta w$", r"$(d\Delta V/dx)/\Delta\rho g$", r"$\Delta$ deflection [m]"
    else:
        wtopo, supp = m.topography(), m.grad_x(m.V()) / DRHOG
        lw_, ls_, ylt = r"topography $w$", r"$(dV/dx)/\Delta\rho g$", "deflection [m]"
    axt.plot(m.xkm[wsel], wtopo[wsel], "C4", lw=4.2, label=lw_)                     # fat solid so the dashed overlay reads
    axt.plot(m.xkm[wsel], supp[wsel], "k", lw=2.6, ls=(0, (7, 4)), label=ls_)
    axt.axhline(0, color="0.7", lw=0.6)
    axt.invert_yaxis()   # z-down: deflection is positive-DOWN, so read 0 at top → positive depth at the bottom (mirrors the warped plate)
    verticals(axt, top_labels=True)
    axt.set_xlim(x0km, Lkm); axt.set_xticklabels([])
    axt.set_ylabel(ylt, fontsize=13); axt.tick_params(labelsize=11)
    axt.text(0.008, 0.85, "(a)", transform=axt.transAxes, fontsize=15, fontweight="bold")
    axt.legend(fontsize=11, loc="lower right", frameon=False, ncol=2)
    fig.add_subplot(gs[0, 1]).axis("off")

    for r, (key, label, cmap) in enumerate(PANELS):
        ax = fig.add_subplot(gs[r + 1, 0])
        ax.imshow(imgs[r], extent=[x0km, Lkm, y0km, y1km], aspect="auto", origin="upper")
        if key != "diff":                   # skip reference lines on the σxx−σzz panel: keep the focus on the principal stresses
            verticals(ax)
        if False:                           # phase-shift annotation commented out (may not make the manuscript)
            tr = blended_transform_factory(ax.transData, ax.transAxes)
            xa, xb = lines["moment_max"], lines["shear_max"]
            ax.annotate("", xy=(xb, 0.13), xytext=(xa, 0.13), xycoords=tr,
                        arrowprops=dict(arrowstyle="<->", color="k", lw=1.5))
            ax.text(0.5 * (xa + xb), 0.18, r"$\lambda/8$" + "\nphase shift", transform=tr, ha="center",
                    va="bottom", fontsize=10, linespacing=0.9,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
            ax.text(0.985, 0.88, r"$\lambda=2\pi\alpha\approx%d$ km" % round(lam / 1e3), transform=ax.transAxes,
                    ha="right", va="top", fontsize=9.5, color="0.2",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
        if key == "rho" and not is_diff:    # signed equivalent density (dipole-moment) centroid depth on the map
            ax.plot(m.xkm[wsel_c], yc_warp[wsel_c] / 1e3, color="black", lw=3.0,
                    label=r"$\hat{\rho}$ centroid",
                    path_effects=[pe.withStroke(linewidth=5.5, foreground="white")])  # halo: visible on dark blue
            jmid = int(np.argmin(np.abs(m.z - m.H / 2)))                          # mid-plate reference layer (z = h/2)
            ymid = -(m.z[jmid] + WARP_FACTOR * f["uy"][:, jmid]) * VEXAG          # warped height of mid-plate (follows the deflection)
            xmid = (m.x + f["ux"][:, jmid]) / 1e3                                  # its displayed x (true horizontal displacement)
            ax.plot(xmid[wsel_c], ymid[wsel_c] / 1e3, color="white", lw=1.6, ls=(0, (5, 3)),
                    label="mid-plate ($h/2$)", zorder=6,
                    path_effects=[pe.withStroke(linewidth=3.0, foreground="black")])  # white dashes, black halo — reads over the centroid
            ax.legend(loc="lower right", fontsize=11, frameon=False, ncol=2)
        ax.set_yticks([]); ax.set_xlim(x0km, Lkm)
        ax.tick_params(labelsize=11)
        ax.set_xticklabels([])
        ax.text(0.008, 0.88, f"({'bcdef'[r]})", transform=ax.transAxes, fontsize=15, fontweight="bold")
        cax = fig.add_subplot(gs[r + 1, 1])
        cb = fig.colorbar(ScalarMappable(Normalize(-clims[r], clims[r]), cmap), cax=cax)
        cb.set_label((r"$\Delta$ " + label) if is_diff else label, fontsize=12); cax.tick_params(labelsize=11)

    # panel (e): resultants, ALL direct stress-field integrals — no analytic/flexure input, no ad-hoc
    # offsets.  N_D=∫(σxx−σzz)dz, the PURE bending moment M=∫σxx(z−h/2)dz, and V=∫σxz dz vanish naturally in
    # the undeflected far field (massless: σ→0), so they are plotted raw.  GPE carries the lithostat, so ΔGPE
    # is referenced.  M is the σxx moment (NOT of σxx−σzz) so dM/dx = V exactly (σzz would break it).
    # MANUSCRIPT CONVENTION: Δ(·) ≡ (·)(x) − (·)(x_T), trench-referenced; ΔGPE* = −Δσ̄zz is zero at the trench and rises
    # to the pull (+2.54 TN/m) at x_I; identity ΔN_D = ΔGPE* (same sign), so N_D (dash-dot) rides ON ΔGPE*.  M is the PURE bending moment ∫σxx(z−z_mid) ⇒ dM/dx (dashed) = V EXACTLY.
    axd = fig.add_subplot(gs[nmap + 1, 0])
    # panel (e): resultants — the ONE definition.  V, ΔN_D, M, ΔGPE* are Cauchy integrals on the DEFORMED
    # column (deformed_line), trench-referenced; NOT the reference-grid integrate_z.  The identity
    # ΔN_D = ΔGPE* means N_D (dash-dot) rides ON ΔGPE*; dM/dx (dashed) = V (pure σxx bending moment).
    xr, Varr, ND, Mn, dGPE = resultants(m)
    xpull = lines["shear_max"]; pull = trench_pull(m)[0]
    axd.plot(xr, Varr / 1e12, "C0", lw=4.0, label=r"$V$")                                      # fat solid
    axd.plot(xr, np.gradient(Mn, xr * 1e3) / 1e12, "k", lw=2.2, ls=(0, (7, 4)), label=r"$dM/dx$")  # ≈ V
    axd.plot(xr, dGPE / 1e12, "C3", lw=4.0, label=r"$\Delta\mathrm{GPE}^{*}$")                   # ΔGPE* (coloured, fat)
    axd.plot(xr, ND / 1e12, "k", lw=2.6, ls=(0, (9, 4, 1.5, 4)), label=r"$\Delta N_D$")         # dash-dot, rides on ΔGPE*
    # the trench pull is ΔGPE* read at the FIRST ISOSTATIC column (first w=0)
    axd.plot(xpull, pull / 1e12, "o", color="C3", ms=10, mec="k", mew=0.9, zorder=7)
    axd.annotate(f"pull = {pull/1e12:.2f} TN m$^{{-1}}$", (xpull, pull / 1e12),
                 textcoords="offset points", xytext=(6, -14), fontsize=10, color="C3", fontweight="bold")
    axd.axhline(0, color="0.7", lw=0.6)
    axd.set_ylabel(r"$V,\ \Delta\mathrm{GPE}^{*},\ \Delta N_D$  [TN m$^{-1}$]", fontsize=12)
    axd.tick_params(labelsize=11)
    axM = axd.twinx()                                  # bending moment on its own axis
    axM.plot(xr, Mn / 1e17, "C2", lw=4.0, label=r"$M$")
    if os.environ.get("HERO_MFIX", "1") == "1":
        axM.set_ylim(-1.75, 0.1)                        # FIXED across the h=60 heroes (contains all M; min ≈ −1.57 for elastic)
    else:
        mlo = Mn.min() / 1e17                           # adaptive (thin plate: M_p ∝ h² is much smaller)
        axM.set_ylim(mlo * 1.18, -mlo * 0.08)
    axM.set_ylabel(r"$M$  [$10^{17}$ N]", color="C2", fontsize=13)
    axM.tick_params(axis="y", labelcolor="C2", labelsize=11)
    verticals(axd)
    axd.set_xlim(x0km, Lkm)
    axd.set_xlabel(r"distance from trench  [km]", fontsize=13)   # x is measured from the trench (x=0)
    axd.text(0.008, 0.88, f"({'bcdef'[nmap]})", transform=axd.transAxes, fontsize=15, fontweight="bold")
    h1, l1 = axd.get_legend_handles_labels(); h2, l2 = axM.get_legend_handles_labels()
    axd.legend(h1 + h2, l1 + l2, ncol=3, fontsize=9.5, loc="lower right", frameon=False, handlelength=3.4)
    fig.add_subplot(gs[nmap + 1, 1]).axis("off")

    if save:
        fig.savefig(out_path, dpi=180)
        print("wrote", out_path)
        rows = [(f"x_{k}", lines[k], "km") for k in REF_KEYS]
        if not is_diff and np.isfinite(zc_depth[wsel_c]).any():
            rows += [("rhohat_centroid_depth_mean_over_window", float(np.nanmean(zc_depth[wsel_c])), "km"),
                     ("rhohat_centroid_depth_min_over_window", float(np.nanmin(zc_depth[wsel_c])), "km"),
                     ("rhohat_centroid_depth_max_over_window", float(np.nanmax(zc_depth[wsel_c])), "km")]
        write_table(os.path.splitext(os.path.basename(out_path))[0], ["quantity", "value", "unit"], rows,
                    script="scripts/render_hero.py", figure=out_path, models=[model_dir] + ([diff_dir] if is_diff else []))
    return fig


if __name__ == "__main__":
    if len(sys.argv) > 1:        # optional: MODEL_DIR  OUT  [DIFF_DIR]
        MODEL_DIR = sys.argv[1]
    if len(sys.argv) > 2:
        OUT = sys.argv[2]
    if len(sys.argv) > 3:
        DIFF_DIR = sys.argv[3]
    build()
