# paper_models.jl — the production driver: every model in the paper is a command of this script.
#
#   suite1       Suite 1: strength models — Tresca baseline, elastic, DD-VM asym, DD-VM sym  → data/suite1_strength
#   nd_sweep     Suite 3: background N_D sweep at V=4                                    → data/suite3_background
#   v_sweep      Suite 2: load V sweep (σ_Y=150)                                       → data/suite2_load
#   thickness    Suite 4: thickness sweep at matched deflection                        → data/suite4_thickness
#   convergence  SI Table S3 convergence of the baseline (own command, not in `all`)   → data/convergence
#   all          suite1 + nd_sweep + v_sweep + thickness
#
#   julia --project=. model/paper_models.jl suite1     # (re)generate one production suite
#   julia --project=. model/paper_models.jl all        # the full publication set
#
# Skip-existing: a model directory that already holds gpe_model.vtu is never touched; a rerun regenerates only what
# is missing (move a finished directory aside to force it).  Every model the driver writes gets a provenance.txt
# stamp (commit, parameters, completion flag).  Production is the MASSLESS platform: no body force, no lithostatic
# prestress, a single Δρg Winkler spring at the top surface, a parabolic end shear at the trench face, an edge
# yield-strength ramp; h = 60 km, σ_Y = 150 MPa (M_p = 1.35e17 N), V = 4 TN/m is the locked baseline.
#
# 2026-09-14: the development-history runs (gravity-on platform, hydrostatic/follower foundation, symmetric springs,
# the exploratory DD-VM variants, the fixed-load thickness members) were removed; regeneration of the baseline, the
# 30-km tuned member and a DD-VM model was verified byte-identical afterwards.

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

# NATIVE CONVENTION: x from the trench (load on the LEFT, x=0); y = DEPTH (y=0 top surface, y=h base).
# Facesets: "left"=trench, "right"=clamp, "bottom"=top surface, "top"=base.
const ρa, ρw, gg, L = 3300.0, 1000.0, 9.81, 1600.0e3
const Δρg = (ρa - ρw) * gg; const Emod = 7.0e10
ramp(bulk) = x -> bulk + 100.0e6 * exp(-x[1] / 10.0e3)         # edge yield-stress ramp at the trench (x=0), λ=10 km

const ELASTIC = 1.0e12     # σ_Y → ∞ ⇒ purely elastic (same setup, for the elastic−plastic difference)
# MASSLESS platform: no body force / prestress / confine (solve_plastic defaults), a single restoring Winkler
# spring at the TOP surface (native faceset "bottom"), stiffness Δρg.  The trench pull lives entirely in the
# spring-induced σzz perturbation of the deflected column (no lithostat to contaminate it).
const sp_massless = [("bottom", Δρg, 0.0)]

"""Uniform massless Tresca (or elastic, σY = ELASTIC) plate of thickness `h` under end shear `V`, with an optional
uniform background in-plane resultant `Nmem` and Young's modulus `E`.  Single top Δρg spring, edge yield-stress ramp,
order 2, no gravity/prestress.  Writes data/<outroot>/<name>."""
function run_tresca_deep(; σY, h, V, nx, nz, name, nsteps = 24, Nmem = 0.0, E = Emod, outroot)
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
        crit = :tresca, springs = sp_massless, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
        confine_x = x -> -Nmem / h, yield = ramp(σY), tract_z = V / h, order = 2,
        nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")   # massless: no gravity/prestress
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-24s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
    write_provenance(outdir; kind = "tresca_deep", stress_frame = "massless", sigma_Y_MPa = σY/1e6, h_km = h/1e3, V_TN = V/1e12,
                     Nmem_TN = Nmem/1e12, E_Pa = E, nx = nx, nz = nz, nsteps = nsteps)
end

"""THICKNESS-SWEEP member (SI Suite 4): uniform massless Tresca at thickness `h`, tuned by V-bisection to a
TARGET trench deflection `w_target` [m] so every plate carries the SAME topography.  Same platform as
run_tresca_deep (single top Δρg spring, edge yield-stress ramp, order 2, no gravity/prestress).  w(V) is
smooth and monotone, so a secant iteration from seed `V0` converges in ~2–3 solves.  Writes to
data/suite4_thickness/<name> and records the converged V
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
                confine_x = x -> 0.0, yield = ramp(σY), tract_z = V / h, order = 2,
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

"""Depth-dependent von Mises (DD-VM): σ_Y(z) = cohesion + μ·ρa·g·depthfac(z) + the trench-edge ramp, with linear
hardening H.  depthfac = z (asymmetric: weak top, strong base) or a tent peaking at mid-plate (`symmetric`).  The
cohesion is chosen so the plastic moment matches the Tresca baseline (M_p = 1.35e17 N).  Massless platform."""
function run_dd_vm_direct(; μ, cohesion, h, V, nx, nz, name, nsteps = 22, H = 7.0e8, symmetric = false, outroot)
    outdir = joinpath(DATA, outroot, name)
    if isdir(outdir) && isfile(joinpath(outdir, "gpe_model.vtu"))
        @printf("%-24s exists — skipping (delete dir to force rerun)\n", name); return
    end
    mkpath(outdir)
    depthfac = symmetric ? (z -> h / 2 - abs(z - h / 2)) : (z -> z)
    ddvmyield = x -> cohesion + μ * ρa * gg * depthfac(x[2]) + 100.0e6 * exp(-x[1] / 10.0e3)   # + trench-edge yield-stress ramp
    ext = symmetric ? cohesion + μ * ρa * gg * h / 2 : cohesion + μ * ρa * gg * h
    @printf("running %-24s DD-VM %-4s μ=%.2f (φ=%.1f°) coh=%.0f %s=%.0f MPa V=%.2f nx=%d nz=%d ...\n",
            name, symmetric ? "sym" : "asym", μ, atand(μ), cohesion / 1e6, symmetric ? "peak" : "base", ext / 1e6, V / 1e12, nx, nz)
    res = solve_plastic(; L = L, h = h, nx = nx, nz = nz, E = Emod, ν = 0.25, σ₀ = cohesion, H = H,
        crit = :mises, springs = sp_massless, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
        confine_x = x -> 0.0, yield = ddvmyield, tract_z = V / h, order = 2,
        nsteps = nsteps, rtol = 1.0e-6, load_face = "left", clamp_face = "right")
    xt, w = topography(res, h; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
    @printf("%-24s trench %5.0f m, yielded %4.1f%%  [OK]\n", name, maximum(w), 100 * ny / length(res.states))
    open(joinpath(outdir, "gpe_topo.csv"), "w") do io
        println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
    end
    export_plastic(res, joinpath(outdir, "gpe_model"))
    write_provenance(outdir; kind = "dd_vm_direct", stress_frame = "massless", mu = μ, cohesion_MPa = cohesion/1e6, symmetric = symmetric,
                     h_km = h/1e3, V_TN = V/1e12, nx = nx, nz = nz, nsteps = nsteps, H = H)
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
        if cmd == "convergence"         # SI Table S3: mesh + load-increment convergence of the baseline (NOT in `all` — the 1200×72 run is expensive)
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
        error("unknown command '$c'. Production commands: suite1 | nd_sweep | v_sweep | thickness | convergence | all.")
    end
end
