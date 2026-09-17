# idealized_beam.jl — Horne's elasto-plastic bending via DISPLACEMENT CONTROL.
#
# A clamped–guided beam (RIGHT edge clamped; LEFT edge a vertical slider pulled down by δ — the native frame,
# load on the left) carries a
# CONSTANT shear V along its length, with the bending moment linear (max at the ends). Driving it by
# the imposed end displacement δ — rather than a load — lets the solve trace *past* the limit load,
# so the outer fibres reach DEEP yield while staying stable, with no foundation and no hardening
# (which is what blurs the shear pinning). There is no transverse distributed load, so σ_zz ≈ 0 —
# Horne's idealized condition — and in the deeply-yielded fibres the shear should pin to ~0, the whole
# shear V being carried by the (parabolic) elastic core.
#
#   * elastic=false → the elasto-plastic case (Tresca, H=0)        → data/idealized_beam/
#   * elastic=true  → benchmark: same geometry, no yield           → data/idealized_beam_elastic/
# A yield-strength ramp protects BOTH ends so yielding develops in the clean span, not the corners.
#
#   julia --project=. model/idealized_beam.jl
include("gpe_plastic.jl")     # J2Plasticity, PState, cell_plastic!, export_plastic, Tresca return (guarded main)

function assemble_vol!(K, g, dh, cv, m, u, states, states_old, yieldfn)
    asm = start_assemble(K, g); nc = getnbasefunctions(cv)
    ke = zeros(nc, nc); ge = zeros(nc)
    nobf = x -> Vec{2}((0.0, 0.0)); nos0 = x -> ZERO_S0₂
    for cell in CellIterator(dh)
        c = cellid(cell); reinit!(cv, cell)
        cell_plastic!(ke, ge, cv, m, u[celldofs(cell)], view(states, :, c), view(states_old, :, c),
                      nobf, nos0, yieldfn, getcoordinates(cell))
        assemble!(asm, celldofs(cell), ke, ge)
    end
    return K, g
end

function solve_beam(; L = 10.0, h = 1.0, nx = 200, nz = 40, E = 60.0, ν = 0.25, σY = 1.0, H = 0.0,
                      δ_max = 2.2, nsteps = 44, order = 2, rtol = 1.0e-8,
                      clamp_boost = 8.0, clamp_λ = 0.05, elastic = false)
    m = J2Plasticity(E, ν, σY, H; crit = :tresca)             # perfect plasticity (H=0); displacement control gives stability
    # yield ramp protects BOTH ends (clamp & guided slider) so yielding develops in the clean span.
    yieldfn = elastic ? (x -> 1.0e6 * σY) :
              (x -> σY * (1.0 + clamp_boost * (exp(-x[1] / (clamp_λ * L)) + exp(-(L - x[1]) / (clamp_λ * L)))))
    grid = generate_grid(Quadrilateral, (nx, nz), Vec(0.0, 0.0), Vec(L, h))
    ip = Lagrange{RefQuadrilateral, order}()^2; qo = 2 * order
    cv = CellValues(QuadratureRule{RefQuadrilateral}(qo), ip)
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    ch = ConstraintHandler(dh)     # NATIVE frame: load (slider) on the LEFT, clamp on the RIGHT; z-down ⇒ pull DOWN is +δ
    add!(ch, Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> [0.0, 0.0], [1, 2]))    # clamp
    add!(ch, Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [+δ_max * t], [2]))      # guided slider, u_z ramped DOWN (+y)
    close!(ch)
    K = allocate_matrix(dh); u = zeros(ndofs(dh)); g = zeros(ndofs(dh))
    nqp = getnquadpoints(cv); ncells = getncells(grid)
    states = [PState() for _ in 1:nqp, _ in 1:ncells]; states_old = [PState() for _ in 1:nqp, _ in 1:ncells]
    for step in 1:nsteps
        t = step / nsteps; update!(ch, t); apply!(u, ch); ref = 1.0    # ramp the prescribed end displacement
        for it in 1:40
            assemble_vol!(K, g, dh, cv, m, u, states, states_old, yieldfn)
            it == 1 && (ref = max(norm(g), 1.0))
            apply_zero!(K, g, ch); nr = norm(g[Ferrite.free_dofs(ch)])
            nr < rtol * ref && break
            it == 40 && error("Newton stalled (step $step), ratio $(nr/ref)")
            u .-= K \ g
        end
        states_old .= states
    end
    return (; grid, dh, u, cv, material = m, states, order, prestress = x -> ZERO_S0₂)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (elastic, name) in ((true, "idealized_beam_elastic"), (false, "idealized_beam"))
        res = solve_beam(; elastic = elastic)
        outdir = joinpath(dirname(@__DIR__), "data", name); mkpath(outdir)
        export_plastic(res, joinpath(outdir, "gpe_model"))
        ny = count(s -> s.k > 0, res.states); κmax = maximum(s.k for s in res.states)
        @printf("%-22s yielded %.1f%%, κ_max=%.4f\n", name, 100 * ny / length(res.states), κmax)
    end
end
