# gpe_mvm.jl — finite-strain (total-Lagrangian, St-Venant–Kirchhoff) 2D plane-strain
# elastic plate, on a Winkler foundation, with a DEAD vertical end load.
#
# This is the finite-strain track toward the GPE model (kept separate from the
# small-strain plasticity2d.jl). MILESTONE 1: with body force off, verify it
# reduces to the analytic beam-on-elastic-foundation flexure (same check as
# verify_flexure.jl). Later milestones add gravity, the hydrostatic follower
# foundation, and plasticity.
#
#   julia --project=. gpe_mvm.jl

using Ferrite, Tensors, LinearAlgebra, Printf

# prestress treatment.
#   :dead    = constant first-PK σ₀ — NO geometric tangent. The lithostat is a reference stress that
#              carries σzz but does NOT modify the plate's bending stiffness ⇒ quasi-thin-plate flexure
#              that matches the analytic BOEF, with a parabolic shear-gradient (∂σxz/∂x first moment ≈ h/2).
#              This is the standard lithospheric-flexure treatment and the DEFAULT.
#   :spatial = follower P₀=J·σ₀·F⁻ᵀ — adds the initial-stress (geometric) stiffness. Because σ₀ is largest
#              at the base it softens the base preferentially, distorting the shear toward the base near
#              the trench (first moment runs to ~h, ΔGPE ~2× the scaling). This effect is far larger than
#              the genuine σ₀/E ~ 2% initial-stress correction and is an ARTIFACT of carrying the lithostat
#              as a stiffness-contributing prestress; kept only for studying the initial-stress effect.
const PRESTRESS_MODE = Ref(:dead)

# ---- constitutive: St-Venant–Kirchhoff (small-strain-accurate finite-strain) ----
struct SVK
    λ::Float64
    μ::Float64
end
function strain_energy(C, m::SVK)
    E = (C - one(C)) / 2
    return m.λ / 2 * tr(E)^2 + m.μ * (E ⊡ E)
end
function constitutive_driver(C, m::SVK)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> strain_energy(y, m), C, :all)
    return 2.0 * ∂Ψ∂C, 2.0 * ∂²Ψ∂C²          # S, ∂S/∂C
end

# ---- element residual + tangent (total Lagrangian) ----
# s0fn(x) returns the SPATIAL (gravity-aligned) Cauchy prestress σ₀(x) at material point x — at
# u=0 (F=I) this equals the lithostat, so u=0 stays equilibrium (no self-compression). It is added
# as a co-rotational FOLLOWER load P₀ = J·σ₀·F⁻ᵀ (Piola transform of a spatial Cauchy stress), which
# keeps the lithostat VERTICAL under material rotation. A frozen material 2nd-PK prestress F·S₀
# instead co-rotates, tilting the huge σzz₀ in the high-curvature trench zone and corrupting the
# stress there (density-dependent dbar/σxz blow-up); the spatial form fixes that. (Wu 2004; Kaus & Becker 2010.)
function assemble_cell!(ke, ge, cv, m, ue, cc, bfn, s0fn)
    n = getnbasefunctions(cv)
    fill!(ke, 0); fill!(ge, 0)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        x = spatial_coordinate(cv, qp, cc)
        bf = bfn(x)                                   # body force at this material point
        σ0 = s0fn(x)                                  # SPATIAL Cauchy prestress at this material point
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u
        C = tdot(F)
        S, ∂S∂C = constitutive_driver(C, m)
        I2 = one(F)
        P = F ⋅ S                                     # elastic first PK
        ∂P∂F = otimesu(I2, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I2)
        if !iszero(σ0)
            if PRESTRESS_MODE[] === :dead             # DEAD initial stress: constant first-PK, NO geometric
                P += σ0                               # tangent ⇒ no initial-stress softening (quasi-thin-plate)
            else                                      # spatial follower: P₀ = J·σ₀·F⁻ᵀ (has geometric stiffness)
                ∂P0∂F, P0 = Tensors.gradient(Fv -> det(Fv) * (σ0 ⋅ transpose(inv(Fv))), F, :all)
                P += P0; ∂P∂F += ∂P0∂F
            end
        end
        for i in 1:n
            δui = shape_value(cv, qp, i)
            ∇δui = shape_gradient(cv, qp, i)
            ge[i] += (∇δui ⊡ P - δui ⋅ bf) * dΩ      # internal − body force
            ∇δuiP = ∇δui ⊡ ∂P∂F
            for j in 1:n
                ∇δuj = shape_gradient(cv, qp, j)
                ke[i, j] += (∇δuiP ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
end

const ẑ = Vec{2}((0.0, 1.0))

function assemble_global!(K, g, dh, cv, fv, m, u, bfn, s0fn, springs, p0, pressset, tractfn, tractset)
    assembler = start_assemble(K, g)
    nc = getnbasefunctions(cv)
    ke = zeros(nc, nc); ge = zeros(nc)
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        assemble_cell!(ke, ge, cv, m, u[celldofs(cell)], getcoordinates(cell), bfn, s0fn)
        assemble!(assembler, celldofs(cell), ke, ge)
    end
    nf = getnbasefunctions(fv)
    kef = zeros(nf, nf); gef = zeros(nf)
    # vertical springs: each (set, k, fref) gives a vertical traction (k·u_z − fref).
    # fref is the reference support at u=0 (carry the column weight); k·u_z is the
    # deflection response (buoyant restoring). Same vertical sign on any face.
    for (sset, kf, fref) in springs
        for fc in FacetIterator(dh, sset)
            reinit!(fv, fc); fill!(kef, 0); fill!(gef, 0)
            ue = u[celldofs(fc)]
            for qp in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, qp)
                uz = function_value(fv, qp, ue) ⋅ ẑ
                for i in 1:nf
                    Niz = shape_value(fv, qp, i) ⋅ ẑ
                    gef[i] += (kf * uz - fref) * Niz * dΓ
                    for j in 1:nf
                        kef[i, j] += kf * Niz * (shape_value(fv, qp, j) ⋅ ẑ) * dΓ
                    end
                end
            end
            assemble!(assembler, celldofs(fc), kef, gef)
        end
    end
    # constant-pressure base p0 (dead upward traction on pressset; no tangent)
    if pressset !== nothing
        for fc in FacetIterator(dh, pressset)
            reinit!(fv, fc); fill!(gef, 0)
            for qp in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, qp)
                for i in 1:nf
                    gef[i] -= p0 * (shape_value(fv, qp, i) ⋅ ẑ) * dΓ
                end
            end
            assemble!(assembler, celldofs(fc), zeros(nf, nf), gef)
        end
    end
    # dead end traction tractfn(x) on the right (vertical shear + horizontal confining
    # = σ₀·n that balances the prestress); position-dependent, residual only.
    for fc in FacetIterator(dh, tractset)
        reinit!(fv, fc); fill!(gef, 0)
        cc = getcoordinates(fc)
        for qp in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, qp)
            tract = tractfn(spatial_coordinate(fv, qp, cc))
            for i in 1:nf
                gef[i] -= (shape_value(fv, qp, i) ⋅ tract) * dΓ
            end
        end
        assemble!(assembler, celldofs(fc), zeros(nf, nf), gef)
    end
    return K, g
end

const ZERO_S0 = zero(SymmetricTensor{2, 2})
function solve(; L, h, nx, nz, E, ν, kfound = 0.0, tract_z = 0.0, found_ref = 0.0,
                 spring_face = "bottom", press_face = nothing, p0 = 0.0,
                 springs = nothing,            # vector of (face_name, k, fref); overrides the single-spring path
                 prestress = x -> ZERO_S0,     # initial 2nd-PK stress σ₀(x) (lithostatic ⇒ no self-compression)
                 confine_x = x -> 0.0,         # horizontal confining traction on the right end = σ₀ₓₓ(z)
                 shear_parabolic = false,      # end shear: false = uniform; true = St-Venant parabola (0 at free surfaces)
                 moment = 0.0,                 # end bending moment M [N·m/m]: linear σxx on the right face, ∫σxx(z−h/2)dz = M
                 order = 1,                    # displacement interpolation order: 1=Q1 (bilinear), 2=Q2 (biquadratic — far better stresses & ∂σxz/∂x, no shear locking)
                 left_dofs = [1, 2],           # which left-face dofs to clamp: [1,2]=full clamp; [2]=z only (plate free to expand in x ⇒ σxx→0, no in-plane softening); a single corner node is then pinned in x to kill rigid translation
                 bodyforce = x -> Vec{2}((0.0, 0.0)), nsteps = 4, rtol = 1.0e-8, verbose = false)
    μ = E / (2(1 + ν)); λ = E * ν / ((1 + ν) * (1 - 2ν))
    material = SVK(λ, μ)
    grid = generate_grid(Quadrilateral, (nx, nz), Vec(0.0, 0.0), Vec(L, h))
    ip = Lagrange{RefQuadrilateral, order}()^2
    qorder = 2 * order                        # enough to integrate the order-p stiffness
    cv = CellValues(QuadratureRule{RefQuadrilateral}(qorder), ip)
    fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(qorder), ip)
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> zeros(length(left_dofs)), left_dofs))
    if !(1 in left_dofs)                   # x free on the left ⇒ pin one corner node in x to remove rigid horizontal translation
        corner = first(i for i in 1:getnnodes(grid) if norm(grid.nodes[i].x - Vec(0.0, 0.0)) < 1.0e-6)
        add!(ch, Dirichlet(:u, Set([corner]), (x, t) -> [0.0], [1]))
    end
    close!(ch)
    K = allocate_matrix(dh); u = zeros(ndofs(dh)); g = zeros(ndofs(dh))
    # resolve the spring list: explicit `springs` kwarg, else the single-spring kwargs
    spring_spec = springs === nothing ?
        (kfound != 0 ? [(spring_face, kfound, found_ref)] : Tuple{String,Float64,Float64}[]) :
        [(string(s[1]), Float64(s[2]), Float64(s[3])) for s in springs]
    springsets = [(getfacetset(grid, nm), k, fr) for (nm, k, fr) in spring_spec]
    tractset = getfacetset(grid, "right")
    pressset = press_face === nothing ? nothing : getfacetset(grid, press_face)
    for step in 1:nsteps
        f = step / nsteps                     # ramp load, gravity, prestress, springs, base pressure together
        bfn = x -> f * bodyforce(x); pp = f * p0; s0fn = x -> f * prestress(x)
        # end loading on the right face = classic shear V + moment M (thin-plate trench BC).
        # Shear σxz(z): uniform, or the St-Venant parabola 1.5(V/h)(1−ξ²), ξ=(z−h/2)/(h/2), which
        # is 0 at the free top/bottom surfaces. Moment: linear σxx(z)=12M(z−h/2)/h³ (∫σxx(z−h/2)dz=M).
        hmid = h / 2; Iz = h^3 / 12
        shfn = shear_parabolic ? (z -> 1.5 * tract_z * (1 - ((z - hmid) / hmid)^2)) : (z -> tract_z)
        mfn = z -> moment * (z - hmid) / Iz                       # bending fibre stress from M
        tractfn = x -> Vec{2}((f * (confine_x(x) + mfn(x[2])), f * shfn(x[2])))   # right end: (confine+moment, shear)
        springs_f = [(set, k, f * fr) for (set, k, fr) in springsets]   # ramp the reference support too
        apply!(u, ch)
        ref = 1.0
        for it in 1:25
            assemble_global!(K, g, dh, cv, fv, material, u, bfn, s0fn, springs_f, pp, pressset, tractfn, tractset)
            # characteristic force = full residual incl. reactions (≫0 even when u=0 is the
            # solution, e.g. the prestressed gravity-only state where the free-dof residual ~0)
            it == 1 && (ref = max(norm(g), 1.0))
            apply_zero!(K, g, ch)
            nr = norm(g[Ferrite.free_dofs(ch)])
            verbose && @printf("    step %d it %d  nr = %.4e  ratio %.4e\n", step, it, nr, nr / ref)
            nr < rtol * ref && break
            it == 25 && error("Newton did not converge (step $step), residual ratio $(nr/ref)")
            u .-= K \ g
        end
    end
    # prestress returned so stress recovery can add it: the genuine Cauchy stress is
    # σ = (1/J)·F·(S(C)+S₀)·Fᵀ. For massless runs prestress≡ZERO_S0 (no change).
    return (; grid, dh, u, cv, material, prestress)
end

# Cauchy σ_zz at the quadrature points, with their current z-coordinate
function szz_profile(res)
    cv, m, dh, u = res.cv, res.material, res.dh, res.u
    zs = Float64[]; szz = Float64[]
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            ∇u = function_gradient(cv, qp, ue)
            Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m)
            σ = (1 / det(Fd)) * Fd ⋅ S ⋅ Fd'      # Cauchy stress
            push!(zs, spatial_coordinate(cv, qp, cc)[2]); push!(szz, σ[2, 2])
        end
    end
    return zs, szz
end

# band-averaged column resultants: N = ∫σxx dz and GPE = ∫σzz dz (Cauchy), over xlo<x<xhi
function band_NG(res, xlo, xhi)
    cv, m, dh, u = res.cv, res.material, res.dh, res.u
    Nsum = 0.0; Gsum = 0.0
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, cc)
            (xlo ≤ x[1] ≤ xhi) || continue
            dΩ = getdetJdV(cv, qp)
            ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m)
            σ = (1 / det(Fd)) * Fd ⋅ S ⋅ Fd'
            Nsum += σ[1, 1] * dΩ; Gsum += σ[2, 2] * dΩ
        end
    end
    w = xhi - xlo
    return Nsum / w, Gsum / w          # ∫σxx dz and ∫σzz dz, averaged over the band
end

# σzz(z) profiles of the left and right columns, binned by CURRENT height (so the
# topographic step shows; σzz = 0 above a column's surface ≡ zero-density air).
function szz_columns(res, xL, xR, hh; nbins = 30)
    cv, m, dh, u = res.cv, res.material, res.dh, res.u
    sL = zeros(nbins); cL = zeros(Int, nbins); sR = zeros(nbins); cR = zeros(Int, nbins)
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            xr = spatial_coordinate(cv, qp, cc)
            zc = xr[2]                                          # reference height (both bases at z=0)
            ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m)
            σzz = ((1 / det(Fd)) * Fd ⋅ S ⋅ Fd')[2, 2]
            b = clamp(Int(floor(zc / hh * nbins)) + 1, 1, nbins)
            if xL[1] ≤ xr[1] ≤ xL[2]; sL[b] += σzz; cL[b] += 1
            elseif xR[1] ≤ xr[1] ≤ xR[2]; sR[b] += σzz; cR[b] += 1
            end
        end
    end
    zb = [(i - 0.5) / nbins * hh for i in 1:nbins]
    σL = [cL[i] > 0 ? sL[i] / cL[i] : 0.0 for i in 1:nbins]    # 0 = above surface (air)
    σR = [cR[i] > 0 ? sR[i] / cR[i] : 0.0 for i in 1:nbins]
    return zb, σL, σR
end

# Pressure (−σzz) along x at a fixed CURRENT spatial level z_ref, plus at the
# deflected base. Tests the user's hypothesis: pressure at a reference equipotential
# is ~constant across columns, while the true (deflected) base pressure varies.
# Returns (xc, p_ref, p_base, w_base): x-bin centres, −σzz at z_ref, −σzz at the
# lowest sampled point in the column, and the basal deflection u_z.
function pressure_vs_x(res, z_ref, hh; nxbins = 40)
    cv, m, dh, u = res.cv, res.material, res.dh, res.u
    L = maximum(n.x[1] for n in res.grid.nodes)
    # per x-bin: σzz nearest z_ref (current z), σzz at the deepest point, base u_z
    sref = fill(NaN, nxbins); dref = fill(Inf, nxbins)
    sbot = fill(NaN, nxbins); zbot = fill(Inf, nxbins); wbot = zeros(nxbins)
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            xr = spatial_coordinate(cv, qp, cc)
            uq = function_value(cv, qp, ue)
            zc = xr[2] + uq[2]                                   # CURRENT height
            ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m)
            σzz = ((1 / det(Fd)) * Fd ⋅ S ⋅ Fd')[2, 2]
            b = clamp(Int(floor(xr[1] / L * nxbins)) + 1, 1, nxbins)
            if abs(zc - z_ref) < dref[b]; dref[b] = abs(zc - z_ref); sref[b] = σzz; end
            if zc < zbot[b]; zbot[b] = zc; sbot[b] = σzz; wbot[b] = uq[2]; end
        end
    end
    xc = [(i - 0.5) / nxbins * L for i in 1:nxbins]
    return xc, -sref, -sbot, wbot
end

# Top-surface check: per x-bin, the surface deflection w_top (u_z at the highest
# point) and σzz there. The top spring should make σzz_surf = −ρ_w·g·(−w_top) =
# the SEAWATER hydrostatic load on a sunk surface (less than a rock load ρ·g·w).
function top_surface(res, ρ, gg, h; nxbins = 40)
    cv, m, dh, u = res.cv, res.material, res.dh, res.u
    L = maximum(n.x[1] for n in res.grid.nodes)
    szz = fill(NaN, nxbins); ztop = fill(-Inf, nxbins); wtop = zeros(nxbins); gap = zeros(nxbins)
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            xr = spatial_coordinate(cv, qp, cc); uq = function_value(cv, qp, ue)
            zc = xr[2] + uq[2]
            b = clamp(Int(floor(xr[1] / L * nxbins)) + 1, 1, nxbins)
            if zc > ztop[b]
                ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
                S, _ = constitutive_driver(tdot(Fd), m)
                ztop[b] = zc; szz[b] = ((1 / det(Fd)) * Fd ⋅ S ⋅ Fd')[2, 2]
                wtop[b] = uq[2]; gap[b] = h - xr[2]      # reference distance from this qp up to the surface
            end
        end
    end
    xc = [(i - 0.5) / nxbins * L for i in 1:nxbins]
    σsurf = szz .+ ρ * gg .* gap            # continue σzz up through rock to the true surface
    return xc, wtop, σsurf
end

# Top-surface topography: u_z(x) of the nodes on the reference top (z = h), sorted by x.
# Uses evaluate_at_grid_nodes (robust) — the global dof index is NOT 2i per grid node.
function topography(res, h; y_surf = h)
    # free-surface deflection w = u_y sampled at y = y_surf.  z-up mesh: surface at y=h (default);
    # native z-down mesh: surface at y=0 (pass y_surf=0.0), and w is then positive-DOWN.
    uv = Ferrite.evaluate_at_grid_nodes(res.dh, res.u, :u)
    nodes = [res.grid.nodes[i].x for i in 1:getnnodes(res.grid)]
    top = [i for i in eachindex(nodes) if abs(nodes[i][2] - y_surf) < 1.0e-6]
    x = [nodes[i][1] for i in top]; w = [uv[i][2] for i in top]
    o = sortperm(x)
    return x[o], w[o]
end

# Export Cauchy-stress fields + displacement to a .vtu for PyVista/ParaView.
# Warp the mesh by "u" in PyVista to see the genuine (finite-strain) deflection;
# color by sigma_zz / sigma_xx / sigma_xz / pressure. The lithostatic background
# is huge (~-1.6 GPa), so also export the *perturbation* of σzz from the columnwise
# reference -ρ_ref·g·(h−z), which is the GPE-bearing signal.
function export_vtk(res, fname; ρref = 0.0, gg = 9.81, h = 1.0, order = 1, litho_ρg = 0.0)
    cv, m, dh, u, grid, s0fn = res.cv, res.material, res.dh, res.u, res.grid, res.prestress
    qr = QuadratureRule{RefQuadrilateral}(2 * order)        # must match the cv quadrature
    nqp = getnquadpoints(cv); ncells = getncells(grid)
    σxx = zeros(nqp, ncells); σzz = zeros(nqp, ncells); σxz = zeros(nqp, ncells)
    pres = zeros(nqp, ncells); dσzz = zeros(nqp, ncells)
    # DYNAMIC σzz (Wu 2004, GJI 158:401): σzz_dyn = σzz − σzz_hydrostat(z_PHYSICAL) =
    # σzz + ρref·g·(h − z_phys), with z_phys = z_ref + u_z. Subtracting the lithostat referenced to
    # the PHYSICAL (deformed) depth removes the "rigid-sink" contamination — because the mesh is
    # total-Lagrangian, ∫σzz over the material column follows the column as it sinks ~|w|, accruing
    # a spurious −½ρg|w|h. σzz_dyn zeroes the local hydrostat pointwise, so ∫σzz_dyn is the genuine
    # missing-mass ΔGPE that tracks topography to the trench. This is THE field to integrate for GPE.
    # (Also exported: *_flex = elastic Cauchy stress (1/J)·F·S·Fᵀ, the prestress excluded — the
    # lithostat-vs-flexural split, useful for visualisation but NOT a GPE fix on its own.)
    dynzz = zeros(nqp, ncells)
    fxx = zeros(nqp, ncells); fzz = zeros(nqp, ncells); fxz = zeros(nqp, ncells)
    for cell in CellIterator(dh)
        c = cellid(cell); reinit!(cv, cell)
        ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:nqp
            xr = spatial_coordinate(cv, qp, cc); uq = function_value(cv, qp, ue)
            ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m)
            iJ = 1 / det(Fd)
            σe = iJ * Fd ⋅ S ⋅ Fd'                            # elastic / flexural part
            σ = σe + s0fn(xr)                                 # total = flexural + SPATIAL prestress σ₀(x)
            z = xr[2]; zphys = z + uq[2]                       # reference vs PHYSICAL height
            lith = litho_ρg * (h - z)                          # isotropic lithostat superposed (0 = massless)
            σxx[qp, c] = σ[1, 1] - lith; σzz[qp, c] = σ[2, 2] - lith; σxz[qp, c] = σ[1, 2]
            dynzz[qp, c] = σ[2, 2] + ρref * gg * (h - zphys)  # Wu dynamic σzz (rigid-sink removed)
            fxx[qp, c] = σe[1, 1]; fzz[qp, c] = σe[2, 2]; fxz[qp, c] = σe[1, 2]
            pres[qp, c] = -(σxx[qp, c] + σzz[qp, c]) / 2       # flexural pressure + lithostatic ρg(h−z)
            dσzz[qp, c] = σ[2, 2] - (-ρref * gg * (h - z))   # σzz perturbation in the REFERENCE frame (still rigid-sink-contaminated)
        end
    end
    projector = L2Projector(Lagrange{RefQuadrilateral, order}(), grid)
    VTKGridFile(fname, dh) do vtk
        write_solution(vtk, dh, u)
        write_projection(vtk, projector, project(projector, σxx, qr), "sigma_xx [Pa]")
        write_projection(vtk, projector, project(projector, σzz, qr), "sigma_zz [Pa]")
        write_projection(vtk, projector, project(projector, σxz, qr), "sigma_xz [Pa]")
        write_projection(vtk, projector, project(projector, pres, qr), "pressure [Pa]")
        write_projection(vtk, projector, project(projector, dσzz, qr), "dsigma_zz [Pa]")
        write_projection(vtk, projector, project(projector, dynzz, qr), "sigma_zz_dyn [Pa]")
        write_projection(vtk, projector, project(projector, fxx, qr), "sigma_xx_flex [Pa]")
        write_projection(vtk, projector, project(projector, fzz, qr), "sigma_zz_flex [Pa]")
        write_projection(vtk, projector, project(projector, fxz, qr), "sigma_xz_flex [Pa]")
    end
    return fname
end

# Axial resultant N(x)=∫σxx dz and GPE proxy ∫σzz dz per column (Cauchy), binned in x.
function NG_profile(res, h; nxbins = 80)
    cv, m, dh, u, s0fn = res.cv, res.material, res.dh, res.u, res.prestress
    L = maximum(n.x[1] for n in res.grid.nodes)
    N = zeros(nxbins); G = zeros(nxbins)
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, cc); dΩ = getdetJdV(cv, qp)
            ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m)
            σ = (1 / det(Fd)) * Fd ⋅ S ⋅ Fd' + s0fn(x)       # spatial prestress added as Cauchy
            b = clamp(Int(floor(x[1] / L * nxbins)) + 1, 1, nxbins)
            N[b] += σ[1, 1] * dΩ; G[b] += σ[2, 2] * dΩ
        end
    end
    binw = L / nxbins
    xc = [(i - 0.5) / nxbins * L for i in 1:nxbins]
    return xc, N ./ binw, G ./ binw       # ∫σxx dz and ∫σzz dz per column (N/m)
end

# σzz at the fixed GEOMETRIC reference level z=0 (current coords), across x: the pressure-
# continuation check. With the physical split (water unload at top, asthenosphere restoring at
# base) the pressure at the z=0 equipotential is held ~constant even as the base sinks below it —
# small spread ⇒ pressure continuity holds. (z_level is a CURRENT height, not material.)
function szz_at_base(res; nxbins = 60, z_level = 0.0)
    cv, m, dh, u, s0fn = res.cv, res.material, res.dh, res.u, res.prestress
    L = maximum(n.x[1] for n in res.grid.nodes)
    s = fill(NaN, nxbins); best = fill(Inf, nxbins)
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            xr = spatial_coordinate(cv, qp, cc); uq = function_value(cv, qp, ue)
            zc = xr[2] + uq[2]; d = abs(zc - z_level)        # current height
            b = clamp(Int(floor(xr[1] / L * nxbins)) + 1, 1, nxbins)
            if d < best[b]
                ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
                S, _ = constitutive_driver(tdot(Fd), m); best[b] = d
                s[b] = ((1 / det(Fd)) * Fd ⋅ S ⋅ transpose(Fd) + s0fn(xr))[2, 2]
            end
        end
    end
    return s
end

# Current base height z_base(x) (reference bottom nodes displaced by u_z), sorted by x.
function base_topography(res)
    uv = Ferrite.evaluate_at_grid_nodes(res.dh, res.u, :u)
    nodes = [res.grid.nodes[i].x for i in 1:getnnodes(res.grid)]
    zb0 = minimum(n[2] for n in nodes)
    base = [i for i in eachindex(nodes) if abs(nodes[i][2] - zb0) < 1.0e-6]
    x = [nodes[i][1] for i in base]; zc = [nodes[i][2] + uv[i][2] for i in base]
    o = sortperm(x); return x[o], zc[o]
end

# Vertically integrated resultants down to a fixed CURRENT level z_eq (the synthetic
# equipotential), per x-bin over [xlo,xhi]: N=∫σxx dz, V=∫σxz dz, I=∫σzz dz (the three
# integrated stress components), plus the pressure −σzz sampled AT z_eq (constant-pressure
# check). Only quadrature points whose current height ≥ z_eq are integrated.
function resultants_to_level(res, z_eq; xlo, xhi, nxbins = 60)
    cv, m, dh, u, s0fn = res.cv, res.material, res.dh, res.u, res.prestress
    N = zeros(nxbins); V = zeros(nxbins); Izz = zeros(nxbins)
    peq = fill(NaN, nxbins); dpeq = fill(Inf, nxbins)
    for cell in CellIterator(dh)
        reinit!(cv, cell); ue = u[celldofs(cell)]; cc = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            xr = spatial_coordinate(cv, qp, cc); uq = function_value(cv, qp, ue)
            x = xr[1]; zc = xr[2] + uq[2]
            (xlo ≤ x ≤ xhi) || continue
            ∇u = function_gradient(cv, qp, ue); Fd = one(∇u) + ∇u
            S, _ = constitutive_driver(tdot(Fd), m); σ = (1 / det(Fd)) * Fd ⋅ S ⋅ Fd' + s0fn(xr)
            b = clamp(Int(floor((x - xlo) / (xhi - xlo) * nxbins)) + 1, 1, nxbins)
            if zc ≥ z_eq
                dΩ = getdetJdV(cv, qp)
                N[b] += σ[1, 1] * dΩ; V[b] += σ[1, 2] * dΩ; Izz[b] += σ[2, 2] * dΩ
            end
            if abs(zc - z_eq) < dpeq[b]; dpeq[b] = abs(zc - z_eq); peq[b] = -σ[2, 2]; end
        end
    end
    binw = (xhi - xlo) / nxbins
    xc = [xlo + (i - 0.5) / nxbins * (xhi - xlo) for i in 1:nxbins]
    return xc, N ./ binw, V ./ binw, Izz ./ binw, peq
end

# ===================== Lithosphere flexure model (single, iterated) =====================
# Massless elastic plate on a buoyant Δρg foundation, loaded by an end shear (trench pull).
# CONSTANT density ⇒ the lithostat is identical in every column and cancels in any column
# difference, so the plate carries NO body force (no self-compression, no prestress, no GPa
# background). Springs carry the buoyant restoring: ρ_asth·g at the base, −ρ_w·g at the top
# → net Δρ·g. Outputs for the notebook: gpe_model.vtu, gpe_topo.csv, gpe_NG.csv.
if abspath(PROGRAM_FILE) == @__FILE__
    ρa, ρw, gg, hh, L = 3300.0, 1000.0, 9.81, 50.0e3, 1600.0e3
    Ep = 7.0e10 / (1 - 0.25^2); D = Ep * hh^3 / 12; Δρg = (ρa - ρw) * gg
    α = (4D / Δρg)^(1 / 4); I2 = one(SymmetricTensor{2, 2}); elt_order = 2

    # Run one configuration and write its outputs to out/<name>/ (gpe_model.vtu + CSVs),
    # printing the key physics checks. Two physics tracks via the kwargs:
    #   MASSLESS (default): no body force, no prestress, single Δρg top spring → clean standard
    #     flexure (no initial-stress softening). The lithostat cancels columnwise.
    #   GENUINE GRAVITY: bodyforce=−ρg ẑ + VERTICAL-ONLY prestress σzz₀=−ρg(h−z) (makes u=0 the
    #     exact gravity equilibrium: zero self-compression, σ₀ carries the real lithostatic σzz).
    #     Springs: top Δρg restoring + bottom dead support ρgh. The genuine σzz (≈−1.6 GPa +
    #     flexural perturbation) is recovered directly from the solve — NO analytic superposition.
    # ρref sets the dσzz perturbation reference (−ρref·g·(h−z)). The trench load is ALWAYS prescribed
    # by its RESULTANT, not the traction: trench_V is the end shear V=∫σxz dz [N/m] (the parabolic end
    # traction is scaled so its integral equals trench_V, i.e. tract_z=trench_V/h); trench_M is the end
    # bending moment M=∫σxx(z−h/2)dz [N·m/m]. Default V=−1 TN/m (peak σxz=1.5·V/h=30 MPa — a physical
    # trench shear; the old −40 MPa traction gave a 60 MPa peak and an unphysical end concentration).
    function run_config(name; trench_M = 0.0, trench_V = -1.0e12, ρref = 0.0,
                        bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0,
                        confine_x = x -> 0.0, springs = [("top", Δρg, 0.0)])
        outdir = joinpath(dirname(@__DIR__), "data", name); mkpath(outdir)
        res = solve(; L = L, h = hh, nx = 400, nz = 24, E = 7.0e10, ν = 0.25, order = elt_order,
            springs = springs, bodyforce = bodyforce, prestress = prestress, confine_x = confine_x,
            tract_z = trench_V / hh, shear_parabolic = true, moment = trench_M, nsteps = 6)  # V resultant → traction

        xt, w = topography(res, hh)
        open(joinpath(outdir, "gpe_topo.csv"), "w") do io
            println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
        end
        xg, Nx, Gx = NG_profile(res, hh; nxbins = 80)
        iref = argmin(abs.(xg .- 0.15L)); iload = argmin(abs.(xg .- 0.9L))
        ΔN = Nx .- Nx[iref]; ΔG = Gx .- Gx[iref]
        open(joinpath(outdir, "gpe_NG.csv"), "w") do io
            println(io, "x,N,GPE,dN,dGPE")
            for i in eachindex(xg); @printf(io, "%.6e,%.6e,%.6e,%.6e,%.6e\n", xg[i], Nx[i], Gx[i], ΔN[i], ΔG[i]); end
        end
        x_fb = xt[argmax(w)]; X0 = L - 2 * (L - x_fb)
        xb, zb = base_topography(res); z_eq = maximum(zb[(xb .≥ X0) .& (xb .≤ L)])
        xr, Nr, Vr, Ir, peq = resultants_to_level(res, z_eq; xlo = X0, xhi = L, nxbins = 60)
        open(joinpath(outdir, "gpe_resultants.csv"), "w") do io
            println(io, "x,N,V,GPE,p_eq")
            for i in eachindex(xr); @printf(io, "%.6e,%.6e,%.6e,%.6e,%.6e\n", xr[i], Nr[i], Vr[i], Ir[i], peq[i]); end
        end
        export_vtk(res, joinpath(outdir, "gpe_model"); ρref = ρref, gg = gg, h = hh, order = elt_order, litho_ρg = 0.0)
        # physics checks
        @printf("\n=== '%s'  V=%.2f TN/m  M=%.1e N·m/m  → data/%s/ ===\n", name, trench_V / 1e12, trench_M, name)
        @printf("  topography: forebulge %+.1f m, trench %.1f m\n", maximum(w), minimum(w))
        @printf("  near load:  ΔN = %+.4f TN/m,  ΔGPE = %+.4f TN/m,  ΔN/ΔGPE = %.3f (ΔN should be ~0)\n",
            ΔN[iload] / 1e12, ΔG[iload] / 1e12, ΔN[iload] / ΔG[iload])
        sb = szz_at_base(res); dom = 16:55                      # skip clamp & loaded-edge bins
        @printf("  σzz at z=0:  mean %.3f GPa, spread %.2f MPa (small ⇒ pressure held at reference level)\n",
            sum(sb[dom]) / length(dom) / 1e9, (maximum(sb[dom]) - minimum(sb[dom])) / 1e6)
        return res
    end

    @printf("flexural length α = %.0f km   L = %.0f km = %.1f α\n", α / 1e3, L / 1e3, L / α)
    run_config("massless")                                   # clean flexure, no lithostat
    run_config("moment"; trench_M = 5.0e16)                   # classic end bending moment
    # GENUINE coupled gravity: vertical-only lithostatic prestress + gravity body force. σzz is
    # solved, not pasted on; flexure is the real (initial-stress-softened) response. Foundation =
    # DEAD BASE (constant lithostatic support ρgh, k=0) + TOP carrying ALL the deflection restoring
    # (Δρg). The base has NO density contrast (plate = asthenosphere density), so a deflection-
    # dependent base reaction is unphysical — and a pressure-continuation base spring (k=ρa·g) injects
    # a spurious base reaction that makes ∫σzz (ΔGPE) NON-MONOTONIC with topography near the trench.
    # The dead base puts the deflection signal at the top (where the water/rock contrast is) ⇒ ΔGPE
    # stays monotonic to the trench. (Verified: PC base turns over near the trench; dead base does not.)
    grav = x -> Vec{2}((0.0, -ρa * gg))
    # ISOTROPIC lithostat σ₀ = −ρg(h−z)·I — pure pressure, NO deviator, so it carries no shear/Mises
    # stress (won't pre-yield) and ALL deviatoric stress comes from the load. The in-plane component
    # σ₀ₓₓ=−ρg(h−z) requires a matching confining traction on the trench face (the lithostatic stress
    # the rest of the subducting plate pushes back with). :dead ⇒ no initial-stress softening.
    σ0iso = x -> -ρa * gg * (hh - x[2]) * I2
    confine = x -> -ρa * gg * (hh - x[2])
    run_config("gravity"; ρref = ρa, bodyforce = grav, prestress = σ0iso, confine_x = confine,
        springs = [("top", Δρg, 0.0), ("bottom", 0.0, ρa * gg * hh)])
end
