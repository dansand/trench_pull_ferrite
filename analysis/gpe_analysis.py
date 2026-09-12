"""gpe_analysis.py — structured-grid analysis for the finite_strain GPE models.

The Ferrite VTU is a rectangular (Nx x Nz) node grid with uniform spacing, so every
field is a clean 2-D array and every depth-integrated quantity is a Simpson/trapz
integral over the mesh nodes — no `sample_over_line` interpolation needed.

    from gpe_analysis import Model, reference_lines
    mod = Model("data/suite1_strength/tresca_deep_150_60km_V4")
    Cxx, Czz, Cxz = mod.cauchy_fields()   # Cauchy; the VTU's sigma_* are 2nd-PK —        # (Nx, Nz) numpy array -> pcolormesh/imshow
    V   = mod.V()                           # reference-grid Cauchy ∫σxz dz (Nx,) — for LOCATING columns only
    dG, dN, x_I = trench_pull(mod)          # the paper's numbers: deformed-Cauchy columns (deformed_resultants)
    lines = reference_lines(mod)            # dict of named x-locations [km]
"""
import json, os
import numpy as np
import pyvista as pv
try:
    from scipy.integrate import simpson
    _INT = lambda y, xc: simpson(y, x=xc, axis=1)
except Exception:                                   # scipy optional -> trapz fallback
    _INT = lambda y, xc: np.trapz(y, x=xc, axis=1)


class Model:
    """A loaded gpe_model.vtu on its native rectangular node grid."""

    def __init__(self, model_dir):
        # CONVENTION (single working frame — the solve now writes it NATIVELY, so this reader is a
        # straight pass-through):
        #   x → to the right, ORIGIN AT THE TRENCH (loaded edge, x=0);  z ↓ positive = DEPTH below
        #   the top surface (z=0 top, z=h base);  the load/trench sits on the LEFT (x=0).  Deflection
        #   w = u_z is positive-DOWN.  σ_xz>0 in the bending core (the top-face shear traction,
        #   outward normal −ẑ, is −σ_xz < 0).  See paper_models.jl for the setup.
        self.dir = model_dir
        self.m = pv.read(f"{model_dir}/gpe_model.vtu")
        p = self.m.points
        xr, zr = np.round(p[:, 0], 6), np.round(p[:, 1], 6)   # 6 dp: clean for both large (lithospheric) and small (non-dim beam) coords
        self.x, self.z = np.unique(xr), np.unique(zr)        # node coords [m]: x from trench, z = depth
        self.Nx, self.Nz = self.x.size, self.z.size
        self.L, self.H = self.x.max(), self.z.max()
        self._ix = np.searchsorted(self.x, xr)
        self._iz = np.searchsorted(self.z, zr)
        self.xkm = self.x / 1e3
        self.depth = self.z / 1e3                            # z IS depth below the top surface [km]

    # ---- fields as structured (Nx, Nz) arrays ----
    def array(self, name, comp=None):
        """Scalar field 'name' as a (Nx, Nz) array. comp picks a vector component
        (0/1/2 or 'x'/'y'/'z'); e.g. array('u','z') is the vertical displacement."""
        v = self.m[name]
        if comp is not None:
            v = v[:, {"x": 0, "y": 1, "z": 2}.get(comp, comp)]
        a = np.full((self.Nx, self.Nz), np.nan)
        a[self._ix, self._iz] = v
        return a                                          # stresses invariant under the 180° frame rotation; u already negated at load

    def names(self):
        return list(self.m.array_names)

    # ---- depth integrals over the mesh nodes -> (Nx,) profiles ----
    def integrate_z(self, field2d):
        return _INT(field2d, self.z)

    # ---- MATERIAL-section resultant: the two idealized-beam benchmarks, nothing else in the paper ----------------
    def resultant_material(self, name):
        """∫ S dz of an exported (2nd Piola–Kirchhoff, material-frame) stress field along the REFERENCE section.
        The right measure for the material-section beam benchmarks (scripts/render_benchmark.py,
        render_mp_benchmark.py: Hetényi V, shear parabola, M–κ on the reference thickness) and used by nothing
        else — every production number goes through the deformed-Cauchy extractor (deformed_line /
        deformed_resultants) below."""
        return self.integrate_z(self.array(name))

    # ---- Reference-grid resultants of the CAUCHY stress: for LOCATING columns only ------------------------------
    # σ = J⁻¹F S Fᵀ is formed at the reference nodes and integrated along the REFERENCE vertical (z at fixed X):
    # a hybrid (spatial stress, material line), exact to O(slope), and adequate for finding where V = 0 or w = 0
    # (reference_lines) and for the w_τ = −(dV/dx)/Δρg panels, which is all the paper uses these for.  No quoted
    # column value comes from them: the trench pull and every column value go through deformed_line /
    # deformed_resultants.  The VTU's sigma_* are 2nd Piola–Kirchhoff: never integrate them raw here (the
    # benchmark exception is resultant_material above).
    def N(self):   return self.integrate_z(self.cauchy_fields()[0])                    # ∫σxx dz  (reference grid)
    def V(self):   return self.integrate_z(self.cauchy_fields()[2])                    # ∫σxz dz  (reference grid)
    def GPE(self): return self.integrate_z(self.cauchy_fields()[1])                    # RAW ∫σzz dz (σzz<0)  (reference grid)
    # SIGN CONVENTION (main.tex): the manuscript GPE ≡ −σ̄_zz = −∫σzz dz = −GPE().  The trench-pull identity
    # is ΔN_D = −Δσ̄_zz ≡ ΔGPE, i.e. ΔN_D = ΔGPE (same sign), both = −∫δσzz.  So for the manuscript ΔGPE use
    # −(GPE()−far); it is NEGATIVE at the trench (pressure deficit) and POSITIVE at the outer rise. fx() = N_D.
    def fx(self):  Cxx, Czz, _ = self.cauchy_fields(); return self.integrate_z(Cxx - Czz)      # N_D = ∫(σxx−σzz)dz  (reference grid)
    def M(self):                                                                       # bending moment ∫σxx(z−h/2)dz  (reference grid)
        return self.integrate_z(self.cauchy_fields()[0] * (self.z - self.H / 2)[None, :])  # z = depth, arm about mid-plane

    def topography(self):
        # NB: this is a 2-D (x,y) mesh — the VERTICAL displacement is component 1 (VTK 'y'),
        # already negated at load, so it is u_z (positive-DOWN).  Top surface is depth index 0.
        return self.array("u", 1)[:, 0]                                                # deflection positive-down [m]

    # ==== BLESSED deformed-frame integration (massless platform; see FINDINGS §16) ====
    # The one audited way to compute the resultants: deform the mesh, form the CAUCHY stress, interpolate it
    # onto a vertical line (structured, layer-by-layer), integrate. Reproduces ΔN_D = ΔGPE* to <0.1%.
    def stress_frame(self):
        """Which stress measure the stored sigma_* fields carry:
          'massless' — 2nd-Piola-Kirchhoff with no prestress (σ0=0); the J⁻¹FSFᵀ push-forward to Cauchy is
                       exact. All production suites are massless.
          'gravity'  — 2nd-Piola + lithostatic prestress (S+σ0); the naive push-forward is inconsistent by
                       the prestress-rotation term (~2–5%) — use the exporter's own *_cauchy fields instead.
        Read from provenance.txt (stress_frame=...) when stamped, else from the central data/DATA_MANIFEST.json
        (the shipped models); otherwise inferred from the σzz scale (lithostatic ≫ the load-induced stress of a
        massless run), with a printed notice."""
        if getattr(self, "_frame", None) is None:
            frame = None
            try:
                with open(f"{self.dir}/provenance.txt") as f:
                    for line in f:
                        if line.startswith("stress_frame="):
                            frame = line.strip().split("=", 1)[1]
            except OSError:
                pass
            if frame is None:                                   # no stamp: the central manifest (shipped models)
                frame = self.manifest_entry().get("stress_frame")
            if frame is None:                                   # neither: infer from the σzz magnitude, and say so
                frame = "gravity" if np.nanmedian(np.abs(self.array("sigma_zz [Pa]"))) > 300e6 else "massless"
                print(f"[gpe_analysis] {self.dir}: no provenance.txt and no DATA_MANIFEST.json entry — stress frame "
                      f"'{frame}' inferred from the σzz magnitude (300 MPa heuristic)")
            self._frame = frame
        return self._frame

    def manifest_entry(self):
        """This model's entry in data/DATA_MANIFEST.json (the central provenance record of the shipped models:
        command, parameters, stress frame, file hashes, origin), found by walking up from the model directory;
        {} if there is none."""
        if getattr(self, "_manifest", None) is None:
            self._manifest = {}
            d = os.path.abspath(self.dir)
            for _ in range(4):
                parent = os.path.dirname(d)
                mf = os.path.join(parent, "DATA_MANIFEST.json")
                if os.path.isfile(mf):
                    with open(mf) as f:
                        models = json.load(f).get("models", {})
                    key = os.path.relpath(os.path.abspath(self.dir), parent).replace(os.sep, "/")
                    self._manifest = models.get(key, {})
                    break
                d = parent
        return self._manifest

    def cauchy_fields(self):
        """Cauchy stress σ = J⁻¹·F·S·Fᵀ on the structured grid, F = I + ∂u/∂X.  Valid for the MASSLESS
        platform, where the exported stress S carries no prestress.  Returns (Cxx, Czz, Cxz), each (Nx,Nz)."""
        if getattr(self, "_cauchy", None) is not None:
            return self._cauchy
        if self.stress_frame() == "gravity":
            raise ValueError(
                f"{self.dir}: cauchy_fields() (and trench_pull / deformed_resultants) are only valid for "
                "MASSLESS exports, where the stored sigma_* fields are 2nd-Piola (σ0=0). This is a GRAVITY "
                "export (sigma_* = S + lithostatic prestress); the J⁻¹FSFᵀ push-forward is inconsistent by the "
                "prestress-rotation term (~2–5%). Use the exporter's consistent *_cauchy fields, or re-run "
                "massless. (All production suites are massless.)")
        Xr, Zr = self.x, self.z
        ux, uz = self.array("u", 0), self.array("u", 1)
        Sxx, Szz, Sxz = self.array("sigma_xx [Pa]"), self.array("sigma_zz [Pa]"), self.array("sigma_xz [Pa]")
        F00 = 1 + np.gradient(ux, Xr, axis=0); F01 = np.gradient(ux, Zr, axis=1)
        F10 = np.gradient(uz, Xr, axis=0);     F11 = 1 + np.gradient(uz, Zr, axis=1)
        J = F00 * F11 - F01 * F10
        FSxx = F00 * Sxx + F01 * Sxz; FSxz = F00 * Sxz + F01 * Szz
        FSzx = F10 * Sxx + F11 * Sxz; FSzz = F10 * Sxz + F11 * Szz
        Cxx = (FSxx * F00 + FSxz * F01) / J
        Czz = (FSzx * F10 + FSzz * F11) / J
        Cxz = (FSxx * F10 + FSxz * F11) / J
        self._cauchy = (Cxx, Czz, Cxz)
        return self._cauchy

    def deformed_line(self, x_km):
        """Cauchy stress profiles along the DEFORMED vertical line at absolute x [km], via structured
        layer-by-layer interpolation.  For each material layer, find the reference X at which its deformed x
        equals x₀, read z_def and the Cauchy σ there.  Massless ⇒ zero stress outside the plate, so a layer
        the line does not cut simply does not contribute.  Returns a dict of depth(z_def)-sorted arrays:
        z (deformed depth), z_ref (the reference depth each sample came from — the exact z↦z_def map on this
        column), sxx, szz, sxz (Cauchy), n_d = sxx−szz.  This is the single per-line extractor; both the
        resultants (integrate) and the stress-profile figures (plot) go through it, so stresses are always
        extracted the same way."""
        x0 = x_km * 1e3
        Cxx, Czz, Cxz = self.cauchy_fields()
        Xr, Zr = self.x, self.z
        ux, uz = self.array("u", 0), self.array("u", 1)
        Xdef = Xr[:, None] + ux; Zdef = Zr[None, :] + uz
        zref, z, sxx, szz, sxz = [], [], [], [], []
        for j in range(self.Nz):
            xl = Xdef[:, j]
            if not (xl.min() <= x0 <= xl.max()):              # line misses this layer ⇒ outside the plate
                continue
            Xi = np.interp(x0, xl, Xr)                        # reference X where this layer crosses x₀
            zref.append(Zr[j]); z.append(np.interp(Xi, Xr, Zdef[:, j]))
            sxx.append(np.interp(Xi, Xr, Cxx[:, j]))
            szz.append(np.interp(Xi, Xr, Czz[:, j]))
            sxz.append(np.interp(Xi, Xr, Cxz[:, j]))
        z = np.array(z); o = np.argsort(z)
        g = lambda a: np.array(a)[o]
        sxx, szz = g(sxx), g(szz)
        return {"z": z[o], "z_ref": g(zref), "sxx": sxx, "szz": szz, "sxz": g(sxz), "n_d": sxx - szz}

    def has_field(self, name):
        return name in self.m.array_names

    def deformed_line_field(self, x_km, field_name, comp=None):
        """Interpolate an exported nodal field onto the DEFORMED vertical line at absolute x [km], using the
        SAME structured layer-by-layer crossing as deformed_line (so the sample points coincide exactly with
        deformed_line's z-grid).  Returns the field values z-sorted to match deformed_line(x)['z'].  Used to
        read the Julia FE-exported Cauchy-stress gradient (dsxz_dx) straight onto the line — no differencing."""
        x0 = x_km * 1e3
        Fld = self.array(field_name, comp)
        Xr = self.x
        ux = self.array("u", 0)
        Xdef = Xr[:, None] + ux
        z, val = [], []
        Zdef = self.z[None, :] + self.array("u", 1)
        for j in range(self.Nz):
            xl = Xdef[:, j]
            if not (xl.min() <= x0 <= xl.max()):
                continue
            Xi = np.interp(x0, xl, Xr)
            z.append(np.interp(Xi, Xr, Zdef[:, j]))
            val.append(np.interp(Xi, Xr, Fld[:, j]))
        z = np.array(z); o = np.argsort(z)
        return np.array(val)[o]

    def deformed_resultants(self, x_km):
        """∫σzz and N_D = ∫(σxx−σzz) along the DEFORMED vertical line at absolute x [km] (thin wrapper over
        deformed_line — the one shared extractor).  Returns (∫σzz, N_D) in N/m."""
        L = self.deformed_line(x_km)
        return np.trapz(L["szz"], L["z"]), np.trapz(L["n_d"], L["z"])

    def batch_deformed_resultants(self, x0s_km):
        """Vectorized deformed_resultants for MANY lines at once — the SAME definition (deformed Cauchy
        interpolated onto vertical lines, integrated over deformed z), with the per-layer crossing search
        batched into one np.interp per layer.  ×100+ faster than looping deformed_resultants; use for the
        many-line ΔGPE*(x) sweeps (overlay / gpe_compare / hero panel-e).  INTERIOR lines only (every layer
        must be cut — keep the scalar deformed_line for lines hugging the trench edge).  Certified bit-
        equivalent to the looped extractor (handoff §5c).  Returns (∫σzz [array], N_D [array]) in N/m."""
        x0 = np.asarray(x0s_km, float) * 1e3
        Cxx, Czz, _ = self.cauchy_fields(); Cnd = Cxx - Czz
        ux, uz = self.array("u", 0), self.array("u", 1)
        Xdef = self.x[:, None] + ux; Zdef = self.z[None, :] + uz
        nL = x0.size
        Z = np.empty((nL, self.Nz)); SZ = np.empty((nL, self.Nz)); ND = np.empty((nL, self.Nz))
        for j in range(self.Nz):
            Xi = np.interp(x0, Xdef[:, j], self.x)          # reference X where layer j crosses each line
            Z[:, j] = np.interp(Xi, self.x, Zdef[:, j])
            SZ[:, j] = np.interp(Xi, self.x, Czz[:, j])
            ND[:, j] = np.interp(Xi, self.x, Cnd[:, j])
        return np.trapz(SZ, Z, axis=1), np.trapz(ND, Z, axis=1)

    def grad_x(self, prof):
        return np.gradient(prof, self.x)                                               # d/dx of an (Nx,) profile

    def grad_x_field(self, field2d):
        return np.gradient(field2d, self.x, axis=0)                                     # ∂/∂x of a (Nx,Nz) field

    def column(self, field2d, x_km):
        """Depth profile (values vs self.depth) of a (Nx,Nz) field at the column nearest x_km."""
        return field2d[np.argmin(np.abs(self.xkm - x_km)), :]


    # ==== NOT USED BY THE PAPER — legacy diagnostics, kept for completeness =====================================
    # Nothing in scripts/, START_HERE.ipynb, tests/ or reproduce.sh calls these; the paper's column values come from
    # deformed_line / deformed_resultants / trench_pull above.
    def resultant(self, name, comp=None):
        """∫ field dz on the reference grid, raw — a generic helper.  For stresses use resultant_material (S) or the
        Cauchy methods above; never integrate the VTU's sigma_* raw for a production number."""
        return self.integrate_z(self.array(name, comp))

    def box_fields(self, rho_w=1000.0, rho_m=3300.0, g=9.81, datum=0.0, z_c=None, npts=500):
        """σzz AND σxx on ONE shared equipotential box, per column, evaluated at DEFORMED positions.
        The box runs from `datum` (SEA LEVEL, z=0) down to `z_c` (default = the deepest deformed plate
        base; below that there is NO FE stress, only the isotropic assumption, so extending z_c deeper
        changes nothing).  This is a THIN-PLATE model, so beyond the plate the stress is an ISOTROPIC
        approximation (σxx = σzz):
          • WATER  above the deformed surface w(x):  σzz = −ρ_w g z  (hydrostatic, 0 at sea level; σxx=σzz)
          • PLATE  w → deformed base:  the FE σzz and FE σxx
          • MANTLE below the deformed base:  σzz continues lithostatically (ρ_m); σxx = σzz
        Because the caps are isotropic, once σzz is known everywhere so is σxx.  Returns (zg[npts],
        SZZ[Nx,npts], SXX[Nx,npts])."""
        sxx, szz, _ = self.cauchy_fields(); uz = self.array("u", 1)
        zdef = self.z[None, :] + uz
        if z_c is None:
            z_c = float(np.nanmax(zdef[:, -1]))                                        # deepest deformed base
        zg = np.linspace(datum, z_c, npts)
        SZZ = np.empty((self.Nx, npts)); SXX = np.empty((self.Nx, npts))
        for i in range(self.Nx):
            zd = zdef[i]; w, base = zd[0], zd[-1]
            sz = np.interp(zg, zd, szz[i])                                             # plate σzz (FE)
            sz = np.where(zg < w, -rho_w * g * zg, sz)                                 # hydrostatic water cap (σzz=0 at sea level)
            sz = np.where(zg > base, szz[i, -1] - rho_m * g * (zg - base), sz)         # lithostatic cap below the base
            sx = np.interp(zg, zd, sxx[i])                                             # plate σxx (FE)
            sx = np.where((zg >= w) & (zg <= base), sx, sz)                            # caps are ISOTROPIC: σxx = σzz
            SZZ[i] = sz; SXX[i] = sx
        return zg, SZZ, SXX

    def box_integrals(self, **kw):
        """∫σzz, ∫σxx and N_D = ∫(σxx−σzz) per column over the shared box (see box_fields).
        Returns dict: zg, SZZ, SXX (fields) and Szz, Sxx, Nd (Nx,)."""
        zg, SZZ, SXX = self.box_fields(**kw)
        Szz = np.trapz(SZZ, zg, axis=1); Sxx = np.trapz(SXX, zg, axis=1)
        return {"zg": zg, "SZZ": SZZ, "SXX": SXX, "Szz": Szz, "Sxx": Sxx, "Nd": Sxx - Szz}

    def gpe_lab(self, z_c=None, rho_w=1000.0, rho_m=3300.0, g=9.81, datum=0.0, return_zc=False):
        """∫σzz per column over the shared equipotential box (thin wrapper on box_fields).  Manuscript
        GPE ≡ −this; trench pull ΔGPE*(x) = −(gpe_lab(x) − gpe_lab(trench))."""
        zg, SZZ, _ = self.box_fields(rho_w=rho_w, rho_m=rho_m, g=g, datum=datum, z_c=z_c)
        out = np.trapz(SZZ, zg, axis=1)
        return (out, float(zg[-1])) if return_zc else out


def _first_sign_change_from_left(xkm, f, edge_skip_km):
    """x [km] of the first sign change of f(x) scanning RIGHT from the trench (x = 0),
    skipping the first `edge_skip_km` (the loaded-edge boundary layer). NaN if none."""
    inb = np.where(xkm > edge_skip_km)[0]
    s = np.sign(f)
    for k in inb[:-1]:
        if s[k] != s[k + 1] and s[k] != 0 and s[k + 1] != 0:
            # linear-interpolate the crossing between k and k+1
            x0, x1, f0, f1 = xkm[k], xkm[k + 1], f[k], f[k + 1]
            return x0 + (x1 - x0) * (-f0) / (f1 - f0)
    return np.nan


def reference_lines(mod, edge_skip_km=6.0):
    """Key x-locations [km] for the flexure (x measured FROM the trench at x=0), all anchored on
    the (lithostat-free) shear resultant V(x) = ∫σxz dz and the topography:
      trench       : loaded left boundary (x = 0)
      moment_max   : first V = 0 crossing from the trench — the MAX BENDING MOMENT (dM/dx = V = 0)
      isostatic    : first deflection zero-crossing (w = 0, ⇒ dV/dx = 0) — the FIRST ISOSTATIC column,
                     the reference for the trench pull.  "shear_max" is kept as an alias (same location).
      outer_rise   : the forebulge (uplift extremum; deflection is positive-down, so its minimum)
      L1           : trench → shear_max distance (the flexural length scale)
      2L1          : 2·L1
    """
    xkm = mod.xkm
    V, w = mod.V(), mod.topography()
    x_moment = _first_sign_change_from_left(xkm, V, edge_skip_km)             # V = 0  -> max moment
    # FIRST ISOSTATIC COLUMN: the first deflection zero-crossing (w = 0 ⇒ no foundation load ⇒ isostatic).
    # For a Winkler foundation dV/dx = −Δρg·w, so this is also the first dV/dx = 0 — but taken directly from
    # w (robust) rather than from grad_x(V) (noisy; it mis-placed the crossing for the membrane case).
    x_shear = _first_sign_change_from_left(xkm, w, edge_skip_km)              # first w = 0 -> first isostatic column
    x_or = xkm[np.argmin(w)]                                                  # forebulge (uplift = min of positive-down w)
    L1 = x_shear
    return {
        "trench": 0.0,
        "moment_max": x_moment,
        "isostatic": x_shear,       # first isostatic column (first w=0) — the reference for the trench pull
        "shear_max": x_shear,       # alias kept for existing callers (profiles stations); same location
        "outer_rise": x_or,
        "L1 (trench→shear_max)": L1,
        "2·L1": 2 * L1,
    }


def signed_centroid(z, tau):
    """Signed centre of mass ∫z·τ dz / ∫τ dz of a depth profile — for a τzx,x (equivalent-density) distribution, the
    depth at which it sits.  ALWAYS signed (CLAUDE.md): it divides by the net resultant ∫τ, so for a sign-changing profile it can
    fall outside the profile's support and is ill-conditioned when |∫τ| ≪ ∫|τ| (a near-balanced dipole) — callers
    guard that case; a magnitude-weighted ∫z|τ|/∫|τ| is a different quantity and is never used."""
    return np.trapz(z * tau, z) / np.trapz(tau, z)


def trench_ref_km(mod):
    """Absolute x [km] of the leftmost COMPLETE deformed vertical column at the trench.
    The trench material edge (reference X=0) tilts and shifts to deformed x = u_x(0,z); a vertical line
    left of its rightmost excursion would skip the layers whose edge sits to its right — integrating a
    column that pokes out of the deformed body (under-counts ΔGPE*).  So the trench reference is
    max_z[ u_x(0,z) ] — the closest fully-inside vertical line to the deformed trench edge.  (When the
    edge stays at x≤0, e.g. no membrane or tension, this reduces to ≈0 and matches the old x=0 line.)"""
    Xdef0 = mod.x[0] + mod.array("u", 0)[0, :]     # deformed x of the trench material edge, per depth layer
    return Xdef0.max() / 1e3 + 1e-3                # +1 m so the critical layer is safely inside


def deformed_shear_gradient(mod, x_km, delta_km=4.0, force_branch=None):
    """τ_zx,x = ∂σxz/∂x on the DEFORMED vertical line at x [km].  Returns (L0, τ) where
    L0 = mod.deformed_line(x) and τ is on L0['z'].

    BRANCH B (preferred, used automatically when the model carries the FE-exported field): read the
    Julia-computed FE spatial gradient of the Cauchy shear (`dsxz_dx`, formed from the L2-projected
    stress via the mapped FE interpolant) straight onto the deformed line — no Python differencing.
    BRANCH A (fallback for older models without the field): interpolate the Cauchy shear σxz onto
    vertical lines at x±δ (off the deformed mesh) and central-difference.
    `force_branch` ('A'/'B') overrides the auto choice (for cross-checking)."""
    L0 = mod.deformed_line(x_km)
    zc = L0["z"]
    use_b = mod.has_field("dsxz_dx [Pa/m]") if force_branch is None else (force_branch == "B")
    if use_b:                                                 # FE-exported gradient, interpolated onto the line
        tau = mod.deformed_line_field(x_km, "dsxz_dx [Pa/m]")
        return L0, tau
    Lm = mod.deformed_line(x_km - delta_km)
    Lp = mod.deformed_line(x_km + delta_km)
    havem, havep = Lm["z"].size > 0, Lp["z"].size > 0
    if havem and havep:                                       # centred
        tau = (np.interp(zc, Lp["z"], Lp["sxz"]) - np.interp(zc, Lm["z"], Lm["sxz"])) / (2 * delta_km * 1e3)
    elif havep:                                               # forward (near the trench edge)
        tau = (np.interp(zc, Lp["z"], Lp["sxz"]) - L0["sxz"]) / (delta_km * 1e3)
    else:                                                     # backward (near the clamp)
        tau = (L0["sxz"] - np.interp(zc, Lm["z"], Lm["sxz"])) / (delta_km * 1e3)
    return L0, tau


def trench_pull(mod):
    """Trench pull ΔGPE* = −Δσ̄zz and ΔN_D, between the TRENCH and the FIRST ISOSTATIC column x_iso
    (first w=0), via the deformed-Cauchy pipeline (Model.deformed_resultants) — NOT the far field (which
    would include the forebulge).  The trench line is the leftmost COMPLETE deformed column (trench_ref_km),
    NOT the fixed absolute x=0: under edge tilt/axial compression the deformed trench edge moves to x>0, so
    x=0 would integrate a column poking out of the body and under-count the pull.  Identity ΔN_D = ΔGPE*
    holds to <0.05%.  Returns (ΔGPE* [N/m], ΔN_D [N/m], x_iso [km])."""
    xiso = reference_lines(mod)["isostatic"]
    x_tr = trench_ref_km(mod)
    SzzT, NdT = mod.deformed_resultants(x_tr)
    SzzI, NdI = mod.deformed_resultants(xiso)
    return -(SzzI - SzzT), NdI - NdT, xiso


# ---- elastic-core half-thickness c(x) (moved here from the retired render_core_thinning.py; used by scripts/render_core_profiles.py) ----
def core_thickness(m):
    """Elastic-core half-thickness c(x) [m] of a plastic bending model by three estimators, plus d = H/2.
    Returns (c_y, c_m, c_s, d): (a) c_y from the contiguous not-yielded band about the neutral plane; (b) c_m from
    inverting the Tresca moment M = σY(d² − c²/3) with the local plateau σY; (c) c_s from the differential-stress
    slope k at the neutral plane, c = σY/k (self-correcting to c = d where the column is fully elastic)."""
    z, H = m.z, m.H; d = H / 2.0
    Cxx, Czz, _ = m.cauchy_fields(); diff = Cxx - Czz
    yld = m.array("yielded") > 0.5
    jmid = int(np.argmin(np.abs(z - H / 2)))
    c_y = np.full(m.Nx, np.nan)
    for i in range(m.Nx):
        el = ~yld[i]
        if not el[jmid]:
            c_y[i] = 0.0; continue
        lo = jmid
        while lo - 1 >= 0 and el[lo - 1]:
            lo -= 1
        hi = jmid
        while hi + 1 < m.Nz and el[hi + 1]:
            hi += 1
        c_y[i] = 0.5 * (z[hi] - z[lo])
    M = m.integrate_z(diff * (z - H / 2)[None, :])
    sigY = np.nanmax(np.abs(diff), axis=1)
    yielded_col = yld.any(axis=1)
    c_m = np.full(m.Nx, np.nan)
    ok = yielded_col & (sigY > 1e6)
    arg = 3.0 * (d ** 2 - np.abs(M) / np.where(sigY > 1e6, sigY, np.nan))
    c_m[ok] = np.sqrt(np.clip(arg[ok], 0.0, None))
    c_s = np.full(m.Nx, np.nan)
    for i in range(m.Nx):
        f = diff[i]; sy = np.nanmax(np.abs(f))
        if sy < 1e6:
            continue
        lin = np.abs(f) < 0.7 * sy
        if lin.sum() < 3:
            continue
        k = np.polyfit(z[lin], f[lin], 1)[0]
        if abs(k) > 1e-6:
            c_s[i] = min(sy / abs(k), d)
    return c_y, c_m, c_s, d


# ==== NOT USED BY THE PAPER — legacy diagnostic, kept for completeness ===========================================
def stress_depth_profiles(mod, line_keys=("moment_max", "shear_max", "outer_rise", "2·L1"),
                          x_window=None):
    """Depth profiles of σxx−σzz, σxz, and ∂σxz/∂x at the named reference lines.
    Returns (fig, lines_dict). Uses matplotlib only (no pyvista plotting)."""
    import matplotlib.pyplot as plt
    lines = reference_lines(mod)
    sxx, szz, sxz = mod.cauchy_fields()
    quantities = [(r"$\sigma_{xx}-\sigma_{zz}$ [MPa]", (sxx - szz) / 1e6),
                  (r"$\sigma_{xz}$ [MPa]", sxz / 1e6),
                  (r"$\partial\sigma_{xz}/\partial x$ [kPa/m]", mod.grad_x_field(sxz) / 1e3)]
    fig, axes = plt.subplots(1, len(quantities), figsize=(4.2 * len(quantities), 4.6), sharey=True,
                             constrained_layout=True)
    for ax, (lab, Q) in zip(axes, quantities):
        for key in line_keys:
            xk = lines[key]
            ax.plot(mod.column(Q, xk), mod.depth, lw=1.7, label=f"{key} ({xk:.0f} km)")
        ax.axvline(0, c="0.6", lw=0.7); ax.set_xlabel(lab); ax.grid(alpha=0.25)
    axes[0].set_ylabel("depth [km]"); axes[0].set_ylim(mod.H / 1e3, 0)
    axes[0].legend(fontsize=7.5, loc="best")
    fig.suptitle(f"{mod.dir}: stress depth profiles at the flexure reference lines", fontsize=11)
    return fig, lines
