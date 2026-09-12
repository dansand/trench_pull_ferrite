"""render_core_profiles.py — τ_zx,x depth profiles across the yielded band (Tresca): FE vs analytic.

Shows the ∩ (load) → ∪ (core-narrowing) morph of the shear-gradient profile inside the elastic core.
The analytic core profile is
    τ_zx,x(ζ) = (3V'/4c³)(c²−ζ²)        [load, ∩, ∫=V']
              + (3Vc'/4c⁴)(3ζ²−c²)      [core-narrowing, ∪, ∫=0]
with ζ measured from the neutral plane.  The z²-coefficient is (−3V'/4c³ + 9Vc'/4c⁴); its sign flips
the parabola:  ∪ (edge-peaked) when  −V' + 3Vc'/c > 0,  else ∩ (centre-peaked).
V, V', c, c' are all measured from the numerical model (V=∫σxz dz smooth; c from moment inversion).

    /opt/anaconda3/envs/pyvista-env/bin/python render_core_profiles.py [MODEL_DIR]
"""
import sys; sys.path.insert(0, "analysis")
import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")   # headless only when run as a script; importable in Jupyter without hijacking the backend
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from gpe_analysis import Model, deformed_shear_gradient
from gpe_analysis import core_thickness

MODEL = "data/suite1_strength/tresca_deep_150_60km_V4"   # symmetric-Tresca baseline (M_p=1.35, V=4)
STATIONS = [60, 90, 120, 150, 175]       # km from trench, spanning the yielded band (V=4: M_max 90, V_max 182)
SMOOTH_KM = 10.0
DELTA_KM = 4.0                            # half-stencil for the deformed-frame ∂σxz/∂x (branch-A fallback only)
OUT = "data/reference_figures/core_profiles_deep60.png"


def fe_tau_deformed(m, xk, d_km):
    """FE τ_zx,x = ∂σxz/∂x on the DEFORMED vertical line at x=xk [km].  Delegates to
    deformed_shear_gradient, which uses the Julia FE-exported gradient (dsxz_dx, branch B) when the
    model carries it and falls back to ±δ central-differencing (branch A) otherwise.  Returns the
    central line's dict L0, its deformed depth zc, and τ(zc) [Pa/m]."""
    L0, tau = deformed_shear_gradient(m, xk, delta_km=d_km)
    return L0, L0["z"], tau


def main():
    m = Model(MODEL); z = m.z; H = m.H; dx = m.x[1] - m.x[0]; xkm = m.xkm
    sig = SMOOTH_KM * 1e3 / dx
    V = m.V(); Vs = gaussian_filter1d(V, sig); Vp = gaussian_filter1d(V, sig, order=1) / dx
    c_y, c_m, c_sl, d = core_thickness(m)                       # use the slope-method c (smooth, edge-stable)
    good = np.isfinite(c_sl); ci = np.interp(np.arange(m.Nx), np.where(good)[0], c_sl[good])
    cs = gaussian_filter1d(ci, sig); cp = gaussian_filter1d(ci, sig, order=1) / dx
    sxz = m.array("sigma_xz [Pa]"); txx_old = m.grad_x_field(sxz)      # OLD: material-frame ∂S_xz/∂X (for comparison only)
    Cxx, Czz, _ = m.cauchy_fields(); diff = Cxx - Czz

    def neutral(i):                                   # z where σxx−σzz = 0, nearest mid-depth
        f = diff[i]; jm = int(np.argmin(np.abs(z - H / 2)))
        for dj in range(m.Nz // 2):
            for j in (jm - dj, jm + dj):
                if 0 <= j < m.Nz - 1 and f[j] * f[j + 1] < 0:
                    return z[j] + (z[j + 1] - z[j]) * (-f[j]) / (f[j + 1] - f[j])
        return H / 2

    fig, axes = plt.subplots(1, len(STATIONS), figsize=(2.8 * len(STATIONS), 4.8),
                             sharey=True, constrained_layout=True)
    print(f"{'x[km]':>6} {'shape':>8}  {'new-vs-old %':>12}  {'vs-equilibrium %':>16}   (τ_zx,x extraction check)")
    for ax, xk in zip(axes, STATIONS):
        i = int(np.argmin(np.abs(xkm - xk)))
        zn = neutral(i); ca, Va, Vpa, cpa = cs[i], Vs[i], Vp[i], cp[i]
        zeta = z - zn; core = np.abs(zeta) < ca
        ana = (3 * Vpa / (4 * ca ** 3)) * (ca ** 2 - zeta ** 2) + (3 * Va * cpa / (4 * ca ** 4)) * (3 * zeta ** 2 - ca ** 2)
        ana_full = np.where(core, ana, 0.0)          # analytic: parabola in the core, ZERO in the plastic wings
        # FE τ_zx,x from the DEFORMED-Cauchy pipeline: Cauchy σxz interpolated onto vertical lines at xk±δ,
        # central-differenced; plotted against DEFORMED depth. The reference-frame constructs (neutral plane,
        # core bounds, analytic curve) are mapped z_ref↦z_def by this column's own map so all share one axis.
        L0, zc, tau = fe_tau_deformed(m, xk, DELTA_KM)
        tau = gaussian_filter1d(tau, 1.0)                       # light z-smoothing of the differenced derivative (branch A)
        dmap = lambda zr: np.interp(zr, L0["z_ref"], L0["z"])   # reference depth -> deformed depth, this column
        d = dmap(z)                                             # reference nodes on the deformed axis
        fe = np.interp(d, zc, tau)                              # differenced τ, aligned to the reference layers
        ax.plot(fe / 1e3, d / 1e3, "C0", lw=1.3, ls=(0, (3, 2)), alpha=0.75)               # full FE (dashed through the wings)
        ax.plot(np.where(core, fe, np.nan) / 1e3, d / 1e3, "C0", lw=2.8, label="FE (Cauchy, deformed)")
        ax.plot(ana_full / 1e3, d / 1e3, "C3", lw=2.0, ls=(0, (6, 3)), label="analytic")    # parabola→0 (drop = the jump)
        ax.axvline(0, c="0.7", lw=0.6); ax.axhline(dmap(zn) / 1e3, c="0.8", lw=0.6, ls=":")
        for zb in (zn - ca, zn + ca):
            ax.axhline(dmap(zb) / 1e3, c="0.85", lw=0.6)   # elastic-core boundaries (on the deformed axis)
        coeff = -Vpa + 3 * Va * cpa / ca              # >0 ⇒ ∪ edge-peaked
        ax.set_title(f"x={xk} km   {'∪ edge' if coeff > 0 else '∩ centre'}\n"
                     f"c={ca/1e3:.0f} km, c'={cpa*1e3:+.2f}, V={Va/1e12:+.2f}", fontsize=8.5)
        ax.set_xlabel(r"$\tau_{zx,x}$ [kPa/m]"); ax.grid(alpha=0.25)
        # diagnostics over the elastic core: new(Cauchy,deformed) vs old(2nd-PK,material) and vs equilibrium −∂σzz/∂z
        tau_eq = gaussian_filter1d(-np.gradient(L0["szz"], L0["z"]), 1.0)           # current-config equilibrium (massless), same smoothing
        core_zc = np.abs(zc - dmap(zn)) < ca                                        # elastic core on the central line
        scale = np.max(np.abs(fe[core])) or 1.0
        d_old = 100 * np.max(np.abs(fe[core] - txx_old[i][core])) / scale
        d_eq = 100 * np.max(np.abs(tau[core_zc] - tau_eq[core_zc])) / scale
        print(f"{xk:6.0f} {'∪ edge' if coeff > 0 else '∩ centre':>8}  {d_old:12.1f}  {d_eq:16.1f}")
    axes[0].set_ylabel("deformed depth  [km]"); axes[0].set_ylim(H / 1e3, 0)
    axes[0].legend(fontsize=9, loc="lower left")
    fig.suptitle(r"$\tau_{zx,x}$ across the yielded band — FE vs analytic (load $\cap$ + core-narrowing $\cup$), Tresca",
                 fontsize=11)
    fig.savefig(OUT, dpi=170); print("wrote", OUT)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        MODEL = sys.argv[1]
    if len(sys.argv) > 2:      # optional stations (comma-separated km) for a different flexure length
        STATIONS = [float(s) for s in sys.argv[2].split(",")]
    if len(sys.argv) > 3:
        OUT = sys.argv[3]
    main()
