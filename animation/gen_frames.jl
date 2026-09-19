# animation/gen_frames.jl — STANDALONE, DISPOSABLE.  Generates per-load-step states for the
# trench-pull LOADING animation.  It reuses ../gpe_plastic.jl (the solver) READ-ONLY and does NOT
# modify any core script.  The whole animation/ directory is safe to delete outright.
#
# Approach: rather than hook the solver's internal load loop (which would mean editing the core),
# we solve the reference Tresca model at a sequence of load fractions f·V (f = 0.12 … 1.0).  The
# loading is monotonic and proportional, so the converged state at load f·V reproduces that point
# in the ramp — good enough for an illustrative animation.  Coarse mesh (disposable, not publication).
#
#   julia --project=. animation/gen_frames.jl   (from the repo root; a model run: 24 solves)        # writes animation/frames/step_NN/gpe_model.vtu
#
# Skip-existing: delete animation/frames to force a rebuild.

include(joinpath(@__DIR__, "..", "model", "gpe_plastic.jl"))
using Printf

# --- reference massless Tresca config (copied here so this script depends on NOTHING but the solver) ---
const Lx   = 1600.0e3
const Hpl  = 60.0e3
const Emod = 7.0e10
const Δρg  = (3300.0 - 1000.0) * 9.81
const spring   = [("bottom", Δρg, 0.0)]            # single top Winkler (massless platform)
const σY   = 150.0e6
const Vfin = 4.0e12                                # final end-shear resultant
ramp(bulk) = x -> bulk + 100.0e6 * exp(-x[1] / 10.0e3)   # trench-edge yield-stress ramp

const NX = 400          # coarse mesh — this is a throwaway animation, not a paper figure
const NZ = 24
const NSTEPS = 20
const NFRAMES = 24

const OUTROOT = joinpath(@__DIR__, "frames")
mkpath(OUTROOT)

fracs = collect(range(0.12, 1.0; length = NFRAMES))
for (k, f) in enumerate(fracs)
    V = f * Vfin
    d = joinpath(OUTROOT, @sprintf("step_%02d", k))
    if isfile(joinpath(d, "gpe_model.vtu"))
        @printf("frame %02d exists — skip\n", k); continue
    end
    mkpath(d)
    @printf("frame %02d/%d  V=%.2f TN/m ...\n", k, NFRAMES, V / 1e12)
    try
        res = solve_plastic(; L = Lx, h = Hpl, nx = NX, nz = NZ, E = Emod, ν = 0.25, σ₀ = σY, H = 0.0,
            crit = :tresca, springs = spring, bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂,
            confine_x = x -> 0.0, yield = ramp(σY), tract_z = V / Hpl, order = 2,
            nsteps = NSTEPS, rtol = 1.0e-6, load_face = "left", clamp_face = "right")
        xt, w = topography(res, Hpl; y_surf = 0.0); ny = count(s -> s.k > 0, res.states)
        @printf("   w=%.0f m  yielded %.1f%%\n", maximum(w), 100 * ny / length(res.states))
        export_plastic(res, joinpath(d, "gpe_model"))
        open(joinpath(d, "meta.txt"), "w") do io; @printf(io, "V_TN=%.4f\nframe=%d\n", V / 1e12, k); end
    catch e
        @printf("   frame %02d FAILED: %s\n", k, sprint(showerror, e))
    end
end
@printf("done — frames in %s\n", OUTROOT)
