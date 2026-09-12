# paper_models.jl — AUTHORITATIVE production driver for the publication model set.
#
# Renamed from sandbox_M4.jl (2026-08): this file (not the old gravity paper_models.jl, now
# gravity_equivalence_demo.jl) generates the committed publication suites. Production is the MASSLESS deep
# h=60 km study, locked at σ_Y=150 MPa (M_p=1.35), V=4 TN/m.
#
#   PRODUCTION branches (write DIRECTLY to data/<suite>/, with a provenance stamp):
#     suite1     Suite 1: strength models — Tresca baseline, elastic, DD-VM asym, DD-VM sym  → data/suite1_strength
#     nd_sweep   Suite 3: background N_D sweep at V=4                                    → data/suite3_background
#     v_sweep    Suite 2: load V sweep (σ_Y=150)                                       → data/suite2_load
#     thickness  Suite 4: thickness sweep at matched deflection                        → data/suite4_thickness
#     convergence  SI Table S2 convergence of the baseline (own command, not in `all`)   → data/convergence
#     all        suite1 + nd_sweep + v_sweep + thickness
#
#   julia --project=. paper_models.jl suite1  # (re)generate one production suite
#   julia --project=. paper_models.jl all     # the full publication set
#
# Existing model dirs are SKIPPED (never overwritten): a rerun regenerates only what is MISSING. Each model
# dir carries a provenance.txt (git hash, mesh, material, load, nsteps, complete flag). Everything below the
# "DIAGNOSTIC / EXPLORATORY" banner is dev-history (writes to data/_sandbox), retained but not the paper set;
# it is destined for diagnostics/ in the release layout (see CLAUDE.md "Release / what ships").

include("gpe_plastic.jl")
using Printf, Dates

# --- provenance: a small stamp next to each model so a regeneration proves its code + parameters ---
git_short() = try readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`) catch; "unknown" end
const DATA = joinpath(dirname(@__DIR__), "data")   # all model output lives at <repo>/data (this file sits in model/)
function write_provenance(outdir; kw...)
    open(joinpath(outdir, "provenance.txt"), "w") do io
        println(io, "ferrite_commit=", git_short())
        println(io, "generated=", Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))
        for (k, v) in kw; println(io, k, "=", v); end
        println(io, "complete=true")          # written LAST, after the VTU export ⇒ its presence = a finished model
    end
end

# ---- M4 setup (identical to paper_models.jl) ----
# NATIVE CONVENTION: x from the trench (load on the LEFT, x=0); y = DEPTH (y=0 top surface, y=h base);
# gravity DOWN = +y.  Facesets: "left"=trench, "right"=clamp, "bottom"=top surface, "top"=base.
const I2 = one(SymmetricTensor{2, 2})
const ρa, ρw, gg, hh, L = 3300.0, 1000.0, 9.81, 50.0e3, 1600.0e3
const Δρg = (ρa - ρw) * gg; const Emod = 7.0e10
const grav    = x -> Vec{2}((0.0, +ρa * gg))                   # gravity points down = +y
const σ0iso   = x -> -ρa * gg * x[2] * I2                      # lithostat: compressive, ∝ depth y
const confine = x -> +ρa * gg * x[2]                           # σxx x-traction on the left face (n=−x)
const sp_one  = [("bottom", Δρg, 0.0), ("top", 0.0, -ρa * gg * hh)]   # restoring at top surface (y=0); dead base at y=h (fref=−ρg·h)
# Symmetric twin of sp_one for the GRAVITY-OFF control: same total restoring stiffness (Δρg) split
# evenly over both faces, no dead base (no weight to hold). Symmetric restoring ⇒ antisymmetric δσzz
# ⇒ ΔGPE→0.  Gravity does not change the deviatoric/residual structure (M1≡M2), so it is left off.
const sp_sym  = [("bottom", Δρg / 2, 0.0), ("top", Δρg / 2, 0.0)]
# Gravity-ON symmetric: same restoring split, but keep the dead base so the lithostatic weight (and
# its stabilising confining pressure) is retained ⇒ deep yielding is stable, matching M4. Symmetric
# restoring still gives antisymmetric δσzz. (Gravity-off symmetric cannot yield deeply without collapse.)
const sp_sym_grav = [("bottom", Δρg / 2, 0.0), ("top", Δρg / 2, -ρa * gg * hh)]
# HYDROSTATIC / follower-support foundation (the physically-referenced alternative to the dead base):
#   top surface  → WATER LOAD  t_z = +ρw·g·w  (destabilising; the "missing water" where rock deflects down)
#   base         → HYDROSTATIC support  t_z = −(ρa·g·h + ρa·g·w)  (grows with CURRENT depth h+w; follower)
# Net stiffness = −ρw·g + ρa·g = Δρg — identical to sp_one, so the deflection is unchanged. The difference
# is the σzz REFERENCE: σzz becomes hydrostatic on absolute equipotentials (constant σzz at constant depth
# outside the bending zone). Small-slope form: the kf·u_z term is the follower magnitude dependence; normal
# rotation / facet-area change (2nd order at ~1–2° slopes) are neglected. fref unchanged ⇒ same reference state.
const sp_hydro = [("bottom", -ρw * gg, 0.0), ("top", ρa * gg, -ρa * gg * hh)]
const Vload   = +1.5e12                                        # load pulls DOWN (+y) at the trench ⇒ tract_z=+V/h>0
ramp(bulk) = x -> bulk + 100.0e6 * exp(-x[1] / 10.0e3)         # edge yield-stress ramp at the trench (x=0), λ=10 km

"""Run one M4 model at yield stress σY [Pa] and resolution (nx,nz); write to data/_sandbox/<name>.
`springs`/`gravity` select the foundation and body force (default = the asymmetric, gravity-on M4)."""
function run_M4(; σY, nx, nz, name, nsteps = 22, springs = sp_one, gravity = true)
    outdir = joinpath(DATA, "_sandbox", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-22s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    @printf("running %-18s σY=%.0f MPa  nx=%d nz=%d (order 2, gravity=%s) ...\n", name, σY/1e6, nx, nz, gravity)
    bf = gravity ? grav : (x -> Vec{2}((0.0, 0.0)))
    ps = gravity ? σ0iso : (x -> ZERO_S0₂)
    cf = gravity ? confine : (x -> 0.0)
    res = solve_plastic(; L = L, h = hh, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = σY, H = 0.0,
        crit = :tresca, springs = springs, bodyforce = bf, prestress = ps, confine_x = cf,
        yield = ramp(σY), tract_z = Vload / hh, moment = 0.0, order = 2, nsteps = nsteps, rtol = 1.0e-6,
        load_face = "left", clamp_face = "right")
    xt, w = topography(res, hh; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)   # w positive-DOWN
    @printf("%-22s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
end

const ELASTIC = 1.0e12     # σ_Y → ∞ ⇒ purely elastic (same setup, for the elastic−plastic difference)
# MASSLESS elastic: no body force / prestress / confine (solve_plastic defaults), single restoring Winkler
# at the TOP surface (native faceset "bottom"), stiffness Δρg. The trench pull lives entirely in the
# spring-induced σzz perturbation of the deflected column (no lithostat to contaminate it).
const sp_massless = [("bottom", Δρg, 0.0)]
function run_massless(; nx, nz, name, nsteps = 10)
    outdir = joinpath(DATA, "_sandbox", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-22s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    @printf("running %-18s MASSLESS elastic  nx=%d nz=%d (order 2) ...\n", name, nx, nz)
    res = solve_plastic(; L = L, h = hh, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = ELASTIC, H = 0.0,
        crit = :mises, springs = sp_massless, yield = nothing, tract_z = Vload / hh, moment = 0.0,
        order = 2, nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")   # massless: defaults = no bodyforce/prestress/confine
    xt, w = topography(res, hh; y_surf = 0.0)
    @printf("%-22s trench %5.0f m  [OK]\n", name, maximum(w))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
end

"""STRONGER + THICKER massless Tresca, targeting DEEPER deflection: custom thickness `h`, yield `σY`, load
`V`.  A thicker plate (bending capacity ∝ h³) and a higher yield stress carry a larger load WITHOUT plastic
breakdown, so the deflection can be pushed past the ~1.3 km cap of the h=50 km / σY=83 MPa M4 model.  Same
massless platform (single top Δρg spring, no gravity/prestress/confine), edge yield-stress ramp, order 2."""
function run_tresca_deep(; σY, h, V, nx, nz, name, nsteps = 24, Nmem = 0.0, E = Emod, ramp_edge = true, moment = 0.0, shear_parabolic = true, load_uz = nothing, springs = sp_massless, follower = false, sym_V = nothing, sym_patch = 4.0e3, face_sxx = nothing, face_sxz = nothing, outroot = "_sandbox")
    outdir = joinpath(DATA, outroot, name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-24s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    @printf("running %-20s DEEP Tresca  σY=%.0f MPa  h=%.0f km  V=%.2f TN/m  Nmem=%.1f TN/m  nx=%d nz=%d ...\n",
            name, σY/1e6, h/1e3, V/1e12, Nmem/1e12, nx, nz)
    # in-plane resultant Nmem: uniform σxx = Nmem/h applied as the confine-face traction (t_x = −σxx on the
    # left face), reacted by the clamp ⇒ uniform N_D += Nmem (>0 = tension). Massless ⇒ no lithostatic confine.
    res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = E, ν = 0.25, σ₀ = σY, H = 0.0,
        crit = :tresca, springs = springs, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
        confine_x = x -> -Nmem / h, yield = ramp_edge ? ramp(σY) : nothing,
        tract_z = load_uz === nothing ? V / h : 0.0, moment = moment, order = 2,
        shear_parabolic = shear_parabolic, load_uz = load_uz, follower = follower, sym_V = sym_V, sym_patch = sym_patch,
        face_sxx = face_sxx, face_sxz = face_sxz,
        nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")   # massless: no gravity/prestress
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-24s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
    write_provenance(outdir; kind = "tresca_deep", stress_frame = "massless", sigma_Y_MPa = σY/1e6, h_km = h/1e3, V_TN = V/1e12,
                     Nmem_TN = Nmem/1e12, E_Pa = E, nx = nx, nz = nz, nsteps = nsteps, moment = moment)
end

"""THICKNESS-SWEEP member (SI Suite 4): uniform massless Tresca at thickness `h`, tuned by V-bisection to a
TARGET trench deflection `w_target` [m] so every plate carries the SAME topography.  Same platform as
run_tresca_deep (single top Δρg spring, edge yield-stress ramp, order 2, no gravity/prestress).  w(V) is
smooth and monotone, so a secant iteration from seed `V0` converges in ~2–3 solves.  Writes to
data/suite4_thickness/<name> (a STABLE, committed location — not _sandbox) and records the converged V
in tuned_V.txt so the hand-off is reproducible.  The h=60 member is the locked baseline
(suite1_strength/tresca_deep_150_60km_V4, V=4) and is NOT regenerated here."""
function run_tresca_thickness(; σY, h, w_target, V0, nx, nz, name, tol_m = 8.0, maxit = 6, nsteps = 24)
    outdir = joinpath(DATA, "suite4_thickness", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-22s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    @printf("running %-18s THICKNESS Tresca  σY=%.0f MPa  h=%.0f km  target w=%.0f m  seed V=%.2f  nx=%d nz=%d ...\n",
            name, σY/1e6, h/1e3, w_target, V0/1e12, nx, nz)
    # A thin/strong Tresca plate has a plastic-COLLAPSE load only just above the matching load, and w(V)
    # goes near-vertical there, so we must APPROACH FROM BELOW and back off on a stall — never bracket from
    # above (that jumps straight into collapse).  solve_at returns (w,res) or `nothing` if Newton stalls.
    ncall = Ref(0)
    solve_at = V -> begin
        ncall[] += 1
        try
            res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = σY, H = 0.0,
                crit = :tresca, springs = sp_massless, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
                confine_x = x -> 0.0, yield = ramp(σY), tract_z = V / h, moment = 0.0, order = 2,
                nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")
            xt, w = topography(res, h; y_surf = 0.0)
            (maximum(w), res)
        catch e
            @printf("  V=%.3f → STALL (collapse)\n", V/1e12); nothing
        end
    end
    # low anchor (must converge); climb geometrically to bracket w_target, halving the step on any stall
    Vlo = 0.7V0; r = solve_at(Vlo)
    r === nothing && (@printf("%-22s ABORT: low anchor stalled at V=%.3f\n", name, Vlo/1e12); return)
    wlo, reslo = r; @printf("  anchor V=%.3f → w=%.0f m\n", Vlo/1e12, wlo)
    Vhi, whi, reshi = Vlo, wlo, reslo; step = 1.12
    while whi < w_target && ncall[] < 14
        r = solve_at(Vhi * step)
        if r === nothing                                          # collapsed: shrink the increment, don't give up
            step = 1 + (step - 1) / 2
            step < 1.012 && (@printf("%-22s CANNOT MATCH — collapses at w≈%.0f m < target %.0f (V≈%.3f TN/m)\n",
                                     name, whi, w_target, Vhi/1e12); return)
            continue
        end
        Vlo, wlo, reslo = Vhi, whi, reshi
        Vhi, whi, reshi = Vhi * step, r[1], r[2]; @printf("  climb V=%.3f → w=%.0f m\n", Vhi/1e12, whi)
    end
    # secant between the good bracket points (both converged, wlo < target ≤ whi).  Vc is kept in LOCKSTEP
    # with the exported solution resc — the trial load is only committed to Vc after its solve succeeds — so a
    # mid-loop stall can never record a load (tuned_V.txt) that disagrees with the exported VTU.
    Vc, wc, resc = Vhi, whi, reshi
    for it in 1:maxit
        abs(wc - w_target) < tol_m && break
        Vtry = clamp(Vlo + (Vhi - Vlo) * (w_target - wlo) / (whi - wlo), Vlo, Vhi)
        r = solve_at(Vtry); r === nothing && break                       # stall: keep the last good (Vc,wc,resc)
        Vc, wc, resc = Vtry, r[1], r[2]; @printf("  it%-2d  V=%.3f → w=%.0f m\n", it, Vc/1e12, wc)
        wc < w_target ? (Vlo, wlo, reslo = Vc, wc, resc) : (Vhi, whi, reshi = Vc, wc, resc)
    end
    if abs(wc - w_target) >= tol_m                                       # explicit tolerance gate: never export or
        @printf("%-22s NOT CONVERGED: w=%.0f m vs target %.0f (|Δ|=%.0f m > tol %.0f); NOTHING WRITTEN — widen bracket/increase maxit and rerun\n",
                name, wc, w_target, abs(wc - w_target), tol_m); return   # label a miss as CONVERGED; write nothing
    end
    ny = count(s -> s.k > 0, resc.states)
    @printf("%-22s CONVERGED V=%.3f TN/m  trench %5.0f m (target %.0f)  yielded %4.1f%%  [OK]\n",
            name, Vc/1e12, wc, w_target, 100*ny/length(resc.states))
    xt, w = topography(resc, h; y_surf = 0.0)
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    open(joinpath(outdir, "tuned_V.txt"), "w") do io   # committed record of the converged load
        @printf(io, "h_km=%.0f\nsigma_Y_MPa=%.0f\nV_TN=%.4f\nw_target_m=%.1f\nw_achieved_m=%.1f\nyielded_pct=%.2f\n",
                h/1e3, σY/1e6, Vc/1e12, w_target, wc, 100*ny/length(resc.states))
    end
    export_plastic(resc, joinpath(outdir, "gpe_model"))
    write_provenance(outdir; kind = "tresca_thickness", stress_frame = "massless", sigma_Y_MPa = σY/1e6, h_km = h/1e3, V_TN = Vc/1e12,
                     w_target_m = w_target, w_achieved_m = wc, nx = nx, nz = nz, nsteps = nsteps)
end

"""Fixed-load member of the thickness V-SWEEP (subsidiary SI figure): massless Tresca at (h, V) with NO
deflection matching — just solve and record, so the pull–deflection and pull–load paths can be traced for
each thickness.  Same platform as run_tresca_thickness.  Writes to data/suite4_thickness/vsweep/<name>
and stores the applied load in load.txt.  Skip-existing; a Newton stall (collapse) is caught and reported."""
function run_tresca_thick_fixed(; σY, h, V, nx, nz, name, nsteps = 24)
    outdir = joinpath(DATA, "suite4_thickness", "vsweep", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-26s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    @printf("running %-22s h=%.0f km  V=%.3f TN/m  nx=%d nz=%d ...\n", name, h/1e3, V/1e12, nx, nz)
    res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = σY, H = 0.0,
        crit = :tresca, springs = sp_massless, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
        confine_x = x -> 0.0, yield = ramp(σY), tract_z = V / h, moment = 0.0, order = 2,
        nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-26s trench %5.0f m  yielded %4.1f%%  [OK]\n", name, maximum(w), 100*ny/length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    open(joinpath(outdir, "load.txt"), "w") do io; @printf(io, "h_km=%.0f\nV_TN=%.4f\n", h/1e3, V/1e12); end
    export_plastic(res, joinpath(outdir, "gpe_model"))
end

"""Elastic counterpart of run_M4 (same gravity / one-sided spring / load), no yielding.
`springs` selects the foundation (default sp_one = dead base; sp_hydro = hydrostatic follower base)."""
function run_elastic(; nx, nz, name, nsteps = 10, springs = sp_one)
    outdir = joinpath(DATA, "_sandbox", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-22s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    @printf("running %-18s ELASTIC  nx=%d nz=%d (order 2) ...\n", name, nx, nz)
    res = solve_plastic(; L = L, h = hh, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = ELASTIC, H = 0.0,
        crit = :mises, springs = springs, bodyforce = grav, prestress = σ0iso, confine_x = confine,
        yield = nothing, tract_z = Vload / hh, moment = 0.0, order = 2, nsteps = nsteps, rtol = 1.0e-6,
        load_face = "left", clamp_face = "right")
    xt, w = topography(res, hh; y_surf = 0.0)
    @printf("%-22s trench %5.0f m  [OK]\n", name, maximum(w))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
end

"""Depth-dependent von Mises (DD-VM): von Mises with a yield stress that INCREASES with depth,
σ_Y(z) = cohesion + α·ρg·(h−z) + edge ramp (NO cap), otherwise identical to the gravity Tresca M4 (sp_one,
gravity, edge ramp, 2× grid). NOT a true Drucker–Prager — the yield sees only the load-induced deviator
(no pressure/Lode asymmetry), so at a given depth the extension and compression fibres yield at the same
|σxx−σzz|. α sets the depth gradient (calibrated via ζn); a small cohesion gives the free surface strength."""
function run_dd_vm(; cohesion, ζn, nx, nz, name, h = hh, V = Vload, nsteps = 22, H = 7.0e8, Nmem = 0.0, springs = sp_one, gravity = true, symmetric = false)
    outdir = joinpath(DATA, "_sandbox", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-22s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    α = cohesion * (h - 2ζn) / ((ζn^2 - h^2 / 2) * ρa * gg)
    # depth factor: ASYMMETRIC (default) strength ∝ depth y (weak top, strong base); SYMMETRIC (`symmetric`)
    # is a tent peaking at MID-PLATE (a strength envelope: weak both surfaces, strong core) ⇒ symmetric
    # yielding about the neutral plane, like the uniform Tresca.
    depthfac = symmetric ? (z -> h / 2 - abs(z - h / 2)) : (z -> z)
    ddvmyield = x -> cohesion + α * ρa * gg * depthfac(x[2]) + 100.0e6 * exp(-x[1] / 10.0e3)   # edge ramp at the trench (x=0)
    # in-plane MEMBRANE prestress: a uniform σxx = Nmem/h added to the loaded-end confine traction (t_x = −σxx
    # on the left face, n=−x̂), reacted by the clamp ⇒ uniform differential resultant N_D += Nmem. Nmem>0 = tension,
    # which shifts the bending neutral plane DOWN (an alternative to the depth-strength asymmetry). Lithostatic prestress kept.
    dσxx = Nmem / h
    # massless (gravity=false): no body force / lithostatic prestress; keep only the membrane traction on the
    # confine face. The DEPTH-DEPENDENT STRENGTH (ddvmyield) is analytic in depth and carries over UNCHANGED —
    # it never read the model's pressure, so dropping the lithostat does not touch the yield profile.
    bf = gravity ? grav   : (x -> Vec{2}((0.0, 0.0)))
    ps = gravity ? σ0iso  : (x -> ZERO_S0₂)
    cf = gravity ? (x -> confine(x) - dσxx) : (x -> -dσxx)
    @printf("running %-18s DD-VM cohesion=%.0f MPa ζn=%.0f km α=%.3f nx=%d nz=%d V=%.2f Nmem=%.2f TN/m ...\n",
            name, cohesion / 1e6, ζn / 1e3, α, nx, nz, V / 1e12, Nmem / 1e12)
    res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = cohesion, H = H,
        crit = :mises, springs = springs, bodyforce = bf, prestress = ps, confine_x = cf,
        yield = ddvmyield, tract_z = V / h, moment = 0.0, order = 2, nsteps = nsteps, rtol = 1.0e-6,
        load_face = "left", clamp_face = "right")
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-22s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
end

"""DD-VM with DIRECT friction μ (= tanφ) and depth-AVERAGE strength pinned to σavg (= Tresca σY).
σ_Y(z) = cohesion + μ·ρa·g·depthfac(z); <depthfac> over [0,h] = h/2 (asym ramp) or h/4 (sym tent), so
cohesion = σavg − μ·ρa·g·<depthfac>.  Pinning the average pins the plastic moment M_p = σavg·h²/4 (the μ term
integrates out against the bending arm) ⇒ deflection/pull stay ≈ Tresca; the gradient changes only WHICH fibres
yield (weak top / strong base).  Massless platform."""
function run_dd_vm_fric(; μ, σavg, h, V, nx, nz, name, nsteps = 22, H = 7.0e8, Nmem = 0.0, symmetric = false)
    outdir = joinpath(DATA, "_sandbox", name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-24s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    facavg = symmetric ? h / 4 : h / 2
    cohesion = σavg - μ * ρa * gg * facavg                      # pin depth-average = σavg (= Tresca)
    depthfac = symmetric ? (z -> h / 2 - abs(z - h / 2)) : (z -> z)
    ddvmyield = x -> cohesion + μ * ρa * gg * depthfac(x[2]) + 100.0e6 * exp(-x[1] / 10.0e3)
    surf = cohesion; ext = symmetric ? cohesion + μ * ρa * gg * h / 2 : cohesion + μ * ρa * gg * h
    @printf("running %-24s DD-VM-FRIC %-4s μ=%.3f (φ=%.1f°) coh/surf=%.0f  %s=%.0f  avg=%.0f MPa  nx=%d nz=%d ...\n",
            name, symmetric ? "sym" : "asym", μ, atand(μ), surf / 1e6, symmetric ? "peak" : "base", ext / 1e6, σavg / 1e6, nx, nz)
    res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = σavg, H = H,
        crit = :mises, springs = sp_massless, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
        confine_x = x -> -Nmem / h, yield = ddvmyield, tract_z = V / h, moment = 0.0, order = 2,
        nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-24s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
end

"""DD-VM with DIRECT (μ, cohesion): σ_Y(z) = cohesion + μ·ρa·g·depthfac(z).  Used to hit a TARGET plastic
moment M_p (cohesion precomputed for the target) rather than pinning a reference-point strength."""
function run_dd_vm_direct(; μ, cohesion, h, V, nx, nz, name, nsteps = 22, H = 7.0e8, Nmem = 0.0, symmetric = false, edge_ramp = true, outroot = "_sandbox")
    outdir = joinpath(DATA, outroot, name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-24s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    depthfac = symmetric ? (z -> h / 2 - abs(z - h / 2)) : (z -> z)
    rampfac = edge_ramp ? 100.0e6 : 0.0                                   # trench-edge yield-stress ramp (diagnostic toggle)
    ddvmyield = x -> cohesion + μ * ρa * gg * depthfac(x[2]) + rampfac * exp(-x[1] / 10.0e3)
    ext = symmetric ? cohesion + μ * ρa * gg * h / 2 : cohesion + μ * ρa * gg * h
    @printf("running %-24s DD-VM %-4s μ=%.2f (φ=%.1f°) coh=%.0f %s=%.0f MPa V=%.2f nx=%d nz=%d ...\n",
            name, symmetric ? "sym" : "asym", μ, atand(μ), cohesion / 1e6, symmetric ? "peak" : "base", ext / 1e6, V / 1e12, nx, nz)
    res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = cohesion, H = H,
        crit = :mises, springs = sp_massless, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
        confine_x = x -> -Nmem / h, yield = ddvmyield, tract_z = V / h, moment = 0.0, order = 2,
        nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-24s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
    write_provenance(outdir; kind = "dd_vm_direct", stress_frame = "massless", mu = μ, cohesion_MPa = cohesion/1e6, symmetric = symmetric,
                     h_km = h/1e3, V_TN = V/1e12, Nmem_TN = Nmem/1e12, nx = nx, nz = nz, nsteps = nsteps, H = H)
end

if abspath(PROGRAM_FILE) == @__FILE__
    PROD = ("suite1", "nd_sweep", "v_sweep", "thickness", "convergence", "all")
    if !isempty(ARGS) && ARGS[1] in PROD
    # ========================= PRODUCTION (publication suites → data/, provenance-stamped) =========================
        valid = join(PROD, " | ")
        for a in ARGS                        # several commands may be given at once; they run left-to-right
            a in PROD || error("unknown production command '$a'. Valid: $valid.")
        end
        nfail = Ref(0)
        for cmd in ARGS
        if cmd in ("suite1", "all")     # Suite 1: strength models — locked baseline σ_Y=150 (M_p=1.35), V=4
            run_tresca_deep(;  σY = 150e6,   h = 60e3, V = 4e12, nx = 800, nz = 48, outroot = "suite1_strength", name = "tresca_deep_150_60km_V4")
            run_tresca_deep(;  σY = ELASTIC, h = 60e3, V = 4e12, nx = 800, nz = 48, outroot = "suite1_strength", name = "elastic_deep_60km_V4")
            run_dd_vm_direct(; μ = 0.15, cohesion = 31.534e6, h = 60e3, V = 4e12, nx = 800, nz = 48, outroot = "suite1_strength", name = "dd_vm_asym_60km_V4")
            run_dd_vm_direct(; μ = 0.35, cohesion = 36.695e6, h = 60e3, V = 4e12, nx = 800, nz = 48, symmetric = true, outroot = "suite1_strength", name = "dd_vm_sym_60km_V4")
        end
        if cmd in ("nd_sweep", "all")   # Suite 3: background N_D sweep at V=4 (±3; N_D=−4 buckles, +4 dropped for symmetry)
            for Nm in (-3.0, -2.0, -1.0, 1.0, 2.0, 3.0)
                try run_tresca_deep(; σY = 150e6, h = 60e3, V = 4e12, nx = 800, nz = 48, Nmem = Nm*1e12, outroot = "suite3_background", name = "tresca_deep_150_60km_V4_mem$(Int(round(Nm)))")
                catch e; nfail[] += 1; @printf("  N_D=%.0f FAILED: %s\n", Nm, sprint(showerror, e)); end
            end
        end
        if cmd in ("v_sweep", "all")    # Suite 2: load V sweep (σ_Y=150); the V=4 point is the suite1 baseline
            for Vt in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5)
                vn = replace(replace(string(Vt), ".0" => ""), "." => "p")
                try run_tresca_deep(; σY = 150e6, h = 60e3, V = Vt*1e12, nx = 800, nz = 48, outroot = "suite2_load", name = "tresca_deep_150_60km_V$(vn)")
                catch e; nfail[] += 1; @printf("  V=%.1f FAILED: %s\n", Vt, sprint(showerror, e)); end
            end
        end
        if cmd in ("thickness", "all")  # Suite 4: thickness sweep at matched trench deflection (h=60 baseline w=3233 m)
            for (hk, V0k) in ((30.0, 2.05), (40.0, 2.70), (50.0, 3.35))
                nm = "tresca_150_$(Int(round(hk)))km"
                try run_tresca_thickness(; σY = 150e6, h = hk*1e3, w_target = 3232.9, V0 = V0k*1e12, nx = 800, nz = 48, name = nm)
                catch e; @printf("  h=%.0f km FAILED: %s\n", hk, sprint(showerror, e)); end
                # the tuner reports non-convergence by RETURNING (no VTU written), not by throwing — so check
                isfile(joinpath(DATA, "suite4_thickness", nm, "gpe_model.vtu")) ||
                    (nfail[] += 1; @printf("  h=%.0f km INCOMPLETE — no model written (tuner did not converge)\n", hk))
            end
        end
        if cmd == "convergence"         # SI Table S2: mesh + load-increment convergence of the baseline (NOT in `all` — the 1200×72 run is expensive)
            for (nxc, nzc) in ((400, 24), (800, 48), (1200, 72))
                run_tresca_deep(; σY = 150e6, h = 60e3, V = 4e12, nx = nxc, nz = nzc, nsteps = 24, outroot = "convergence", name = "bench_$(nxc)x$(nzc)")
            end
            for ns in (12, 48)
                run_tresca_deep(; σY = 150e6, h = 60e3, V = 4e12, nx = 800, nz = 48, nsteps = ns, outroot = "convergence", name = "bench_800x48_ns$(ns)")
            end
        end
        end                                  # for cmd in ARGS
        if nfail[] > 0
            @printf("\n%d model(s) FAILED to solve — suite INCOMPLETE (see FAILED lines above).\n", nfail[]); exit(1)
        end
    else
        c = isempty(ARGS) ? "" : ARGS[1]
        error("unknown command '$c'. Production commands: suite1 | nd_sweep | v_sweep | thickness | convergence | all.\n" *
              "Diagnostic / exploratory commands live in diagnostics/sandbox.jl.")
    end
end
