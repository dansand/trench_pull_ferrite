# gpe_plastic.jl — J2 (von Mises) elasto-plasticity on the massless flexure platform.
#
# Small-strain J2 in the MATERIAL frame: the return map acts on the Green–Lagrange strain
# E and the 2nd-Piola stress S (both rotation-invariant), so the large flexural rotation is
# carried exactly by the total-Lagrangian kinematics and — since |E| is small — S ≈ σ_Cauchy,
# i.e. this is co-rotational J2 to O(strain). Perfect plasticity (H=0) by default; the surface
# Δρg spring distributes the support so no plastic hinge forms. Reuses the elastic machinery
# (springs/traction/diagnostics) from gpe_mvm.jl. Writes out/plastic/.
#
#   julia --project=. gpe_plastic.jl

include("gpe_mvm.jl")     # topography, base_topography, ẑ, Tensors/Ferrite — main block guarded

# ---- J2 plane-strain return map (from sheet_plasticity/plasticity2d.jl) ----
# σ₀ = mean yield; Δσ_lode = σ_ext − σ_comp encodes a Lode(J3)-dependent strength: the yield stress is
# σ_y(θ) = σ₀ + (Δσ_lode/2)·sin 3θ, so σ_ext at triaxial extension (θ=+30°, the bending tension fibre)
# and σ_comp at compression (θ=−30°, the bottom fibre). Δσ_lode=0 ⇒ ordinary von Mises. Flow is kept
# von-Mises-radial (non-associated), which preserves θ ⇒ the return is consistent (σ_y constant during it).
struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}}
    G::T; K::T; σ₀::T; H::T; Δσ_lode::T; crit::Symbol; Dᵉ::S
end
function J2Plasticity(E, ν, σ₀, H; Δσ_lode = 0.0, crit = :mises)
    δ(i, j) = i == j ? 1.0 : 0.0
    G = E / 2(1 + ν); K = E / 3(1 - 2ν)
    f(i, j, k, l) = 2.0G * (0.5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) + ν / (1.0 - 2.0ν) * δ(i, j) * δ(k, l))
    return J2Plasticity(G, K, σ₀, H, Δσ_lode, crit, SymmetricTensor{4, 3}(f))
end
struct PState{T, S <: SecondOrderTensor{3, T}}
    ϵᵖ::S; σ::S; k::T
end
PState() = PState(zero(SymmetricTensor{2, 3}), zero(SymmetricTensor{2, 3}), 0.0)
vMises(σ) = (s = dev(σ); sqrt(1.5 * s ⊡ s))

# ---- in-plane TRESCA (max in-plane shear) return -------------------------------------------------
# Yields when the in-plane Mohr radius reaches σY/2, i.e. q = √((σxx−σzz)² + 4σxz²) = σY, INDEPENDENT
# of the out-of-plane σyy (which stays elastic). So where σxz≈0 this pins σxx−σzz = ±σY (flat fibre
# stress), unlike von Mises which mixes in σyy. Kludge: radial return on the in-plane deviator only;
# the consistent tangent is taken by AD (Tensors.gradient) of this stress function.
function tresca_stress(ϵ::SymmetricTensor{2, 3}, m::J2Plasticity, st::PState, σY)
    σᵗ = m.Dᵉ ⊡ (ϵ - st.ϵᵖ)
    a = σᵗ[1, 1] - σᵗ[3, 3]; b = 2 * σᵗ[1, 3]; q = sqrt(a^2 + b^2 + 1.0e-18)   # q = √((σxx−σzz)²+4σxz²)
    φ = q - (σY + m.H * st.k)
    φ <= 0.0 && return σᵗ
    μ = φ / (2 * m.G + m.H); f = (q - 2 * m.G * μ) / q                          # scale in-plane deviator
    mxz = (σᵗ[1, 1] + σᵗ[3, 3]) / 2
    return SymmetricTensor{2, 3}((mxz + f * a / 2, σᵗ[1, 2], f * σᵗ[1, 3], σᵗ[2, 2], σᵗ[2, 3], mxz - f * a / 2))
end
function tresca_return(ϵ::SymmetricTensor{2, 3}, m::J2Plasticity, st::PState, σY)
    σ = tresca_stress(ϵ, m, st, σY)
    σᵗ = m.Dᵉ ⊡ (ϵ - st.ϵᵖ); a = σᵗ[1, 1] - σᵗ[3, 3]; q = sqrt(a^2 + (2σᵗ[1, 3])^2 + 1.0e-18)
    φ = q - (σY + m.H * st.k)
    φ <= 0.0 && return σᵗ, m.Dᵉ, PState(st.ϵᵖ, σᵗ, st.k)
    D = Tensors.gradient(e -> tresca_stress(e, m, st, σY), ϵ)                   # consistent tangent via AD
    μ = φ / (2 * m.G + m.H)
    return σ, D, PState(st.ϵᵖ + (σᵗ - σ) / (2 * m.G), σ, st.k + μ)              # Δϵᵖ = Dᵉ⁻¹(σᵗ−σ) = (σᵗ−σ)/2G
end

# ---- BENDING / shear-free fibre return (illustration) --------------------------------------------
# Yields on |σxx−σzz| = σY (σxz NOT in the criterion), AND relaxes the shear σxz→0 in the yielded
# material — the yielded fibre is a shear-free slip layer, so all transverse shear is pushed into the
# elastic core. Gives σxz ≈ 0 wherever yielding occurs (the classic plate "yielded fibres carry only
# σxx" picture). NOTE: zero shear strength ⇒ a soft shear tangent; the elastic core must carry V.
function bending_stress(ϵ::SymmetricTensor{2, 3}, m::J2Plasticity, st::PState, σY)
    σᵗ = m.Dᵉ ⊡ (ϵ - st.ϵᵖ)
    a = σᵗ[1, 1] - σᵗ[3, 3]; q = abs(a); φ = q - (σY + m.H * st.k)
    φ <= 0.0 && return σᵗ
    μ = φ / (2 * m.G + m.H); an = sign(a) * (q - 2 * m.G * μ); mxz = (σᵗ[1, 1] + σᵗ[3, 3]) / 2
    return SymmetricTensor{2, 3}((mxz + an / 2, σᵗ[1, 2], 0.0, σᵗ[2, 2], σᵗ[2, 3], mxz - an / 2))  # σxz → 0
end
function bending_return(ϵ::SymmetricTensor{2, 3}, m::J2Plasticity, st::PState, σY)
    σ = bending_stress(ϵ, m, st, σY)
    σᵗ = m.Dᵉ ⊡ (ϵ - st.ϵᵖ); φ = abs(σᵗ[1, 1] - σᵗ[3, 3]) - (σY + m.H * st.k)
    φ <= 0.0 && return σᵗ, m.Dᵉ, PState(st.ϵᵖ, σᵗ, st.k)
    D = Tensors.gradient(e -> bending_stress(e, m, st, σY), ϵ)
    return σ, D, PState(st.ϵᵖ + (σᵗ - σ) / (2 * m.G), σ, st.k + φ / (2 * m.G + m.H))
end

# trial → yield check → radial return → consistent tangent (small strain, 3D)
function stress_tangent(ϵ::SymmetricTensor{2, 3}, m::J2Plasticity, st::PState, σ0_local = m.σ₀)
    m.crit === :tresca && return tresca_return(ϵ, m, st, σ0_local)
    m.crit === :bending && return bending_return(ϵ, m, st, σ0_local)
    G, H = m.G, m.H
    σᵗ = m.Dᵉ ⊡ (ϵ - st.ϵᵖ); sᵗ = dev(σᵗ); σᵗₑ = sqrt(1.5 * sᵗ ⊡ sᵗ)
    sin3θ = σᵗₑ > 1.0e-3 ? clamp(13.5 * det(sᵗ) / σᵗₑ^3, -1.0, 1.0) : 0.0   # Lode angle of the trial deviator
    σy = σ0_local + 0.5 * m.Δσ_lode * sin3θ                                 # local yield stress (σ0_local carries the depth-dependent von Mises strength)
    φ = σᵗₑ - (σy + H * st.k)
    if φ < 0.0
        return σᵗ, m.Dᵉ, PState(st.ϵᵖ, σᵗ, st.k)
    end
    h = H + 3G; μ = φ / h
    s = (1 - 3G * μ / σᵗₑ) * sᵗ; σ = s + vol(σᵗ); σₑ = σy + H * (st.k + μ)
    δ(i, j) = i == j ? 1.0 : 0.0
    Isd(i, j, k, l) = 0.5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) - δ(i, j) * δ(k, l) / 3
    Q(i, j, k, l) = Isd(i, j, k, l) - 1.5 / σₑ^2 * s[i, j] * s[k, l]
    b = (3G * μ / σₑ) / (1 + 3G * μ / σₑ)
    Dt(i, j, k, l) = -2G * b * Q(i, j, k, l) - 9G^2 / (h * σₑ^2) * s[i, j] * s[k, l]
    D = m.Dᵉ + SymmetricTensor{4, 3}(Dt)
    return σ, D, PState(st.ϵᵖ + 1.5 * μ / σₑ * s, σ, st.k + μ)
end
embed3(ε::SymmetricTensor{2, 2}) = SymmetricTensor{2, 3}((ε[1, 1], 0.0, ε[1, 2], 0.0, 0.0, ε[2, 2]))
extract2(σ::SymmetricTensor{2, 3}) = SymmetricTensor{2, 2}((σ[1, 1], σ[1, 3], σ[3, 3]))
extract2(D::SymmetricTensor{4, 3}) = SymmetricTensor{4, 2}((i, j, k, l) -> D[(1, 3)[i], (1, 3)[j], (1, 3)[k], (1, 3)[l]])

# ---- total-Lagrangian plastic cell: P = F·S(E) + σ₀, tangent via the material elastoplastic D ----
# bfn = body force (gravity); s0fn = ISOTROPIC dead lithostatic prestress (added to P, no geometric
# tangent). σ₀ is pure pressure, so it does NOT enter the deviatoric J2 return map (no pre-yield) —
# the return map sees only the load-induced strain; the lithostat just rides along as σzz background.
function cell_plastic!(ke, ge, cv, m, ue, st, st_old, bfn, s0fn, yieldfn, cc)
    n = getnbasefunctions(cv); fill!(ke, 0); fill!(ge, 0)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp); x = spatial_coordinate(cv, qp, cc)
        bf = bfn(x); σ0 = s0fn(x); σY = yieldfn(x)            # σY: depth-dependent von Mises (DD-VM) yield strength
        ∇u = function_gradient(cv, qp, ue); F = one(∇u) + ∇u
        E = symmetric((tdot(F) - one(F)) / 2)                 # Green–Lagrange (2D, symmetric)
        σ3, D3, st[qp] = stress_tangent(embed3(E), m, st_old[qp], σY)
        S = extract2(σ3); DSE = extract2(D3)                  # 2nd-PK and ∂S/∂E
        I2 = one(F); P = F ⋅ S + σ0                           # + dead isotropic prestress
        ∂P∂F = otimesu(I2, S) + F ⋅ DSE ⊡ otimesu(F', I2)     # ∂S/∂C = DSE/2 ⇒ 2·(…) = DSE
        for i in 1:n
            δui = shape_value(cv, qp, i); ∇δi = shape_gradient(cv, qp, i)
            ge[i] += (∇δi ⊡ P - δui ⋅ bf) * dΩ
            ∇δiP = ∇δi ⊡ ∂P∂F
            for j in 1:n
                ke[i, j] += (∇δiP ⊡ shape_gradient(cv, qp, j)) * dΩ
            end
        end
    end
end

function assemble_plastic!(K, g, dh, cv, fv, m, u, states, states_old, springs, tractfn, tractset, bfn, s0fn, yieldfn; follower = false)
    asm = start_assemble(K, g); nc = getnbasefunctions(cv)
    ke = zeros(nc, nc); ge = zeros(nc)
    for cell in CellIterator(dh)
        c = cellid(cell); reinit!(cv, cell)
        cell_plastic!(ke, ge, cv, m, u[celldofs(cell)], view(states, :, c), view(states_old, :, c), bfn, s0fn, yieldfn, getcoordinates(cell))
        assemble!(asm, celldofs(cell), ke, ge)
    end
    nf = getnbasefunctions(fv); kef = zeros(nf, nf); gef = zeros(nf)
    for (sset, kf, fref) in springs                         # foundation springs
        for fc in FacetIterator(dh, sset)
            reinit!(fv, fc); fill!(kef, 0); fill!(gef, 0); ue = u[celldofs(fc)]
            for qp in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, qp); uz = function_value(fv, qp, ue) ⋅ ẑ
                # direction·area factor `a`: FOLLOWER ⇒ pressure normal to the DEFORMED surface (Nanson a = J F⁻ᵀ N,
                # no spurious tangential traction); otherwise the vertical Winkler dead load (a = ẑ).  Magnitude
                # kf·uz−fref (buoyancy restoring ∝ vertical sinking) is unchanged; only the direction/area map differ.
                local a::Vec{2}
                if follower
                    ∇u = function_gradient(fv, qp, ue); F = one(∇u) + ∇u
                    n = getnormal(fv, qp)
                    a = det(F) * (transpose(inv(F)) ⋅ n)
                    a = (n ⋅ ẑ < 0) ? -a : a          # orient to the vertical (+ẑ) convention so flat ⇒ a = ẑ
                else
                    a = ẑ
                end
                for i in 1:nf
                    di = shape_value(fv, qp, i) ⋅ a
                    gef[i] += (kf * uz - fref) * di * dΓ
                    for j in 1:nf; kef[i, j] += kf * di * (shape_value(fv, qp, j) ⋅ ẑ) * dΓ; end   # approx tangent (pressure-magnitude part)
                end
            end
            assemble!(asm, celldofs(fc), kef, gef)
        end
    end
    for fc in FacetIterator(dh, tractset)                   # end shear (+moment) traction
        reinit!(fv, fc); fill!(gef, 0); cc = getcoordinates(fc)
        for qp in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, qp); tr = tractfn(spatial_coordinate(fv, qp, cc))
            for i in 1:nf; gef[i] -= (shape_value(fv, qp, i) ⋅ tr) * dΓ; end
        end
        assemble!(asm, celldofs(fc), zeros(nf, nf), gef)
    end
    return K, g
end

const ZERO_S0₂ = zero(SymmetricTensor{2, 2})
function solve_plastic(; L, h, nx, nz, E, ν, σ₀, H, springs, tract_z, moment = 0.0, Δσ_lode = 0.0,
                         bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂, confine_x = x -> 0.0,
                         yield = nothing, crit = :mises, shear_parabolic = true, order = 2, nsteps = 10, rtol = 1.0e-8,
                         load_face = "right", clamp_face = "left", load_uz = nothing, follower = false,
                         sym_V = nothing, sym_patch = 4.0e3, face_sxx = nothing, face_sxz = nothing)
    m = J2Plasticity(E, ν, σ₀, H; Δσ_lode = Δσ_lode, crit = crit)
    yieldfn = yield === nothing ? (x -> σ₀) : yield        # depth-dependent von Mises: σ_Y(x) = σ₀ + α·ρg·(h−z); default = constant σ₀
    grid = generate_grid(Quadrilateral, (nx, nz), Vec(0.0, 0.0), Vec(L, h))
    ip = Lagrange{RefQuadrilateral, order}()^2; qo = 2 * order
    cv = CellValues(QuadratureRule{RefQuadrilateral}(qo), ip)
    fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(qo), ip)
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    ch = ConstraintHandler(dh); add!(ch, Dirichlet(:u, getfacetset(grid, clamp_face), (x, t) -> [0.0, 0.0], [1, 2]))
    # DISPLACEMENT-CONTROLLED trench: prescribe the load-face vertical displacement (ramped), u_x free ⇒ shear is the
    # REACTION, no imposed shear shape (bend-to-shape test).  When set, run with tract_z=0 so no traction is applied.
    load_uz !== nothing && add!(ch, Dirichlet(:u, getfacetset(grid, load_face), (x, t) -> [load_uz * t], [2]))
    # SYMMETRY-PLANE control: trench face → symmetry plane (u_x = 0); the shear/σxx face traction is replaced,
    # so causes 1–2 (prescribed-profile mismatch + pointwise σxx=0) cannot exist.  Load enters as a central
    # downward pressure over a short top patch (below).
    sym_V !== nothing && add!(ch, Dirichlet(:u, getfacetset(grid, load_face), (x, t) -> [0.0], [1]))
    close!(ch)
    K = allocate_matrix(dh); u = zeros(ndofs(dh)); g = zeros(ndofs(dh))
    nqp = getnquadpoints(cv); ncells = getncells(grid)
    states = [PState() for _ in 1:nqp, _ in 1:ncells]; states_old = [PState() for _ in 1:nqp, _ in 1:ncells]
    springsets = [(getfacetset(grid, nm), kf, fr) for (nm, kf, fr) in springs]
    tractset = getfacetset(grid, sym_V === nothing ? load_face : "top"); hmid = h / 2; Iz = h^3 / 12
    for step in 1:nsteps
        f = step / nsteps
        bfn = x -> f * bodyforce(x); s0fn = x -> f * prestress(x)
        shfn = shear_parabolic ? (z -> 1.5 * tract_z * (1 - ((z - hmid) / hmid)^2)) : (z -> tract_z)
        tractfn =
            face_sxz !== nothing ?                                                                       # NATURAL-PROFILE reapplication (causes 1-2)
                (x -> Vec{2}((f * (face_sxx === nothing ? 0.0 : face_sxx(x[2])), f * face_sxz(x[2])))) :
            sym_V === nothing ?
                (x -> Vec{2}((f * (confine_x(x) + moment * (x[2] - hmid) / Iz), f * shfn(x[2])))) :       # trench: confine + moment, shear
                (x -> Vec{2}((0.0, x[1] <= sym_patch ? f * (sym_V / sym_patch) : 0.0)))                   # symmetry: central line load as top pressure
        springs_f = [(set, kf, f * fr) for (set, kf, fr) in springsets]
        update!(ch, f); apply!(u, ch); ref = 1.0             # ramp any displacement-controlled load with f=step/nsteps
        for it in 1:40
            assemble_plastic!(K, g, dh, cv, fv, m, u, states, states_old, springs_f, tractfn, tractset, bfn, s0fn, yieldfn; follower = follower)
            it == 1 && (ref = max(norm(g), 1.0))
            apply_zero!(K, g, ch); nr = norm(g[Ferrite.free_dofs(ch)])
            nr < rtol * ref && break
            it == 40 && error("Newton stalled (step $step), ratio $(nr/ref)")
            u .-= K \ g
        end
        states_old .= states                                # commit history after the step converges
    end
    return (; grid, dh, u, cv, material = m, states, order, prestress)
end

# export Cauchy(≈S) fields + von Mises + accumulated plastic strain from the committed states.
# σxx/σzz are the TOTAL stress (return-map deviatoric/elastic part + the isotropic lithostatic
# prestress σ₀); von Mises and plastic_strain are unchanged by σ₀ (isotropic ⇒ no deviator).
function export_plastic(res, fname)
    cv, dh, u, grid, s0fn = res.cv, res.dh, res.u, res.grid, res.prestress; st = res.states
    qr = QuadratureRule{RefQuadrilateral}(2 * res.order); nqp = getnquadpoints(cv); nc = getncells(grid)
    sxx = zeros(nqp, nc); szz = zeros(nqp, nc); sxz = zeros(nqp, nc); pres = zeros(nqp, nc)
    vm = zeros(nqp, nc); κ = zeros(nqp, nc); yld = zeros(nqp, nc)
    sxzC = zeros(nqp, nc); szzC = zeros(nqp, nc)                                # true Cauchy σxz, σzz (for their gradients)
    for cell in CellIterator(dh)
        c = cellid(cell); reinit!(cv, cell); cc = getcoordinates(cell); ue = u[celldofs(cell)]
        for qp in 1:nqp
            x = spatial_coordinate(cv, qp, cc); σ0 = s0fn(x); S = extract2(st[qp, c].σ)
            σ2 = S + σ0                                                         # exported (2nd-PK + lithostat), as before
            sxx[qp, c] = σ2[1, 1]; szz[qp, c] = σ2[2, 2]; sxz[qp, c] = σ2[1, 2]
            pres[qp, c] = -(σ2[1, 1] + σ2[2, 2]) / 2; vm[qp, c] = vMises(st[qp, c].σ); κ[qp, c] = st[qp, c].k
            yld[qp, c] = st[qp, c].k > 0.0 ? 1.0 : 0.0                          # binary: has this point yielded?
            ∇u = function_gradient(cv, qp, ue); F = one(∇u) + ∇u               # consistent Cauchy = J⁻¹(F·S + σ₀)·Fᵀ
            σc = ((F ⋅ S + σ0) ⋅ F') / det(F)                                   # (massless: σ₀=0 ⇒ J⁻¹F·S·Fᵀ)
            sxzC[qp, c] = σc[1, 2]; szzC[qp, c] = σc[2, 2]
        end
    end
    ip = Lagrange{RefQuadrilateral, res.order}()
    proj = L2Projector(ip, grid)
    sxzC_nodal = project(proj, sxzC, qr); szzC_nodal = project(proj, szzC, qr)  # Cauchy σxz, σzz on nodes (proj.dh order)
    # SPATIAL gradients of the Cauchy stress, from the FE interpolant of the projected field pushed to the spatial
    # frame: ∂σ/∂x_i = [F⁻ᵀ · ∇_X σ]_i.  τ_zx,x = ∂σxz/∂x; ∂σzz/∂z is the equilibrium partner (massless: sum=0).
    # Replaces the Python finite-difference ("branch B").
    cvs = CellValues(qr, ip); dsxzdx = zeros(nqp, nc); dszzdz = zeros(nqp, nc)
    for cell in CellIterator(dh)
        c = cellid(cell); reinit!(cv, cell); reinit!(cvs, cell); ue = u[celldofs(cell)]
        sxze = sxzC_nodal[celldofs(proj.dh, c)]; szze = szzC_nodal[celldofs(proj.dh, c)]
        for qp in 1:nqp
            ∇u = function_gradient(cv, qp, ue); Fi = inv(one(∇u) + ∇u)'
            dsxzdx[qp, c] = (Fi ⋅ function_gradient(cvs, qp, sxze))[1]          # ∂σxz/∂x (horizontal) = τ_zx,x
            dszzdz[qp, c] = (Fi ⋅ function_gradient(cvs, qp, szze))[2]          # ∂σzz/∂z (vertical)
        end
    end
    VTKGridFile(fname, dh) do vtk
        write_solution(vtk, dh, u)
        write_projection(vtk, proj, project(proj, sxx, qr), "sigma_xx [Pa]")
        write_projection(vtk, proj, project(proj, szz, qr), "sigma_zz [Pa]")
        write_projection(vtk, proj, project(proj, sxz, qr), "sigma_xz [Pa]")
        write_projection(vtk, proj, sxzC_nodal, "sigma_xz_cauchy [Pa]")        # true Cauchy σxz
        write_projection(vtk, proj, project(proj, dsxzdx, qr), "dsxz_dx [Pa/m]")  # FE spatial gradient of Cauchy σxz
        write_projection(vtk, proj, project(proj, dszzdz, qr), "dszz_dz [Pa/m]")  # FE spatial gradient of Cauchy σzz (equilibrium partner)
        write_projection(vtk, proj, project(proj, pres, qr), "pressure [Pa]")
        write_projection(vtk, proj, project(proj, vm, qr), "von Mises [Pa]")
        write_projection(vtk, proj, project(proj, κ, qr), "plastic_strain")
        write_projection(vtk, proj, project(proj, yld, qr), "yielded")          # 0/1 ⇒ contour at 0.5 for the yield front
    end
    return fname
end

if abspath(PROGRAM_FILE) == @__FILE__
    ρa, ρw, gg, hh, L = 3300.0, 1000.0, 9.81, 50.0e3, 1600.0e3; Δρg = (ρa - ρw) * gg
    Emod = 7.0e10; trench_V = -40.0e6      # 40 MPa trench load (max elastic von Mises 151 MPa)

    # Run one plastic config → out/<name>/ (gpe_model.vtu + gpe_topo.csv) and print checks.
    # σ₀ sets how much yields; H regularizes the localization (band width) — with H=0 the band
    # collapses to the mesh (pathological); a small H makes it mesh-objective. The surface spring
    # owns global stability either way.
    function run_plastic(name; σ₀, H, nsteps = 12)
        outdir = joinpath(@__DIR__, "out", name); mkpath(outdir)
        res = solve_plastic(; L = L, h = hh, nx = 400, nz = 24, E = Emod, ν = 0.25, σ₀ = σ₀, H = H,
            springs = [("top", Δρg, 0.0)], tract_z = trench_V, moment = 0.0, order = 2, nsteps = nsteps)
        xt, w = topography(res, hh)
        open(joinpath(outdir, "gpe_topo.csv"), "w") do io
            println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
        end
        export_plastic(res, joinpath(outdir, "gpe_model"))
        κmax = maximum(s.k for s in res.states); ny = count(s -> s.k > 0, res.states)
        @printf("\n=== '%s'  σ₀=%.0f MPa  H=%.1e  → out/%s/ ===\n", name, σ₀ / 1e6, H, name)
        @printf("  forebulge %+.1f m, trench %.1f m;  κ_max=%.4f, yielded %.1f%%\n",
            maximum(w), minimum(w), κmax, 100 * ny / length(res.states))
        return res
    end

    run_plastic("plastic"; σ₀ = 100.0e6, H = 0.0, nsteps = 12)           # clean, mesh-objective baseline
    run_plastic("plastic_hard"; σ₀ = 80.0e6, H = 7.0e8, nsteps = 15)     # more yielding, regularized (H=E/100)

    # ---- GRAVITY + von Mises ---- isotropic dead lithostat + body force + dead base + trench pressure BC.
    # Von Mises is pressure-independent, so the lithostat does NOT change where it yields (it reproduces the
    # massless yield pattern); the value is (i) proving gravity+plasticity assemble & converge, (ii) the
    # absolute σzz/pressure now rides along for a later pressure- or Lode-dependent criterion. σ₀=50 MPa
    # gives ~0.7× the peak elastic von Mises (~72 MPa) at V=−1 TN/m ⇒ outer-fibre yielding near the trench.
    I2 = one(SymmetricTensor{2, 2})
    grav = x -> Vec{2}((0.0, -ρa * gg))
    σ0iso = x -> -ρa * gg * (hh - x[2]) * I2          # ISOTROPIC pressure ⇒ zero deviator ⇒ no pre-yield
    confine = x -> -ρa * gg * (hh - x[2])             # σ₀ₓₓ confining traction on the trench face
    function run_grav_plastic(name; σ₀, H, V = -1.0e12, M = 0.0, nsteps = 12, yield = nothing, rtol = 1.0e-8)
        outdir = joinpath(@__DIR__, "out", name); mkpath(outdir)
        res = solve_plastic(; L = L, h = hh, nx = 400, nz = 24, E = Emod, ν = 0.25, σ₀ = σ₀, H = H,
            springs = [("top", Δρg, 0.0), ("bottom", 0.0, ρa * gg * hh)],
            bodyforce = grav, prestress = σ0iso, confine_x = confine, yield = yield,
            tract_z = V / hh, moment = M, order = 2, nsteps = nsteps, rtol = rtol)
        xt, w = topography(res, hh)
        open(joinpath(outdir, "gpe_topo.csv"), "w") do io
            println(io, "x,w"); for i in eachindex(xt); @printf(io, "%.6e,%.6e\n", xt[i], w[i]); end
        end
        export_plastic(res, joinpath(outdir, "gpe_model"))
        κmax = maximum(s.k for s in res.states); ny = count(s -> s.k > 0, res.states)
        @printf("\n=== '%s'  σ₀=%.0f MPa  H=%.1e  V=%.2f TN/m  → out/%s/ ===\n", name, σ₀ / 1e6, H, V / 1e12, name)
        @printf("  forebulge %+.1f m, trench %.1f m;  κ_max=%.4f, yielded %.1f%%\n",
            maximum(w), minimum(w), κmax, 100 * ny / length(res.states))
        return res
    end
    run_grav_plastic("gravity_plastic"; σ₀ = 50.0e6, H = 7.0e8, V = -1.0e12, nsteps = 12)

    # ---- GRAVITY + DRUCKER–PRAGER (depth-dependent yield) ---- pressure-dependent strength via the
    # lithostat: σ_Y(ζ) = c + α·ρg·ζ (ζ = hh−z = depth below top). c = COHESION (strength at the free
    # surface ζ=0); α·ρg·ζ = the frictional, pressure-dependent part. Stronger at depth ⇒ the compression
    # (bottom) fibre is stronger than the extension (top) fibre ⇒ the elastic core / neutral plane is pushed
    # DOWN. α set so the fully-plastic plane sits at ζn: solving b·ζn²+2a·ζn=a·h+b·h²/2 (a=c, b=α·ρg) ⇒
    # α = c·(h−2ζn)/((ζn²−h²/2)·ρg). α scales with c, so ζn is invariant to the cohesion magnitude.
    # STABILITY: V=−1.8 over-yielded the end section through-thickness (topography/resultants went anomalous
    # near the trench — a limit-load collapse). Raising the cohesion (c: 30→40 MPa) lifts the whole strength
    # profile and dialling V back (−1.8→−1.2) keeps the section from going fully plastic. ζn stays 30 km.
    let c = 40.0e6, ζn = 30.0e3
        α = c * (hh - 2ζn) / ((ζn^2 - hh^2 / 2) * ρa * gg)
        run_grav_plastic("drucker_prager"; σ₀ = c, H = 7.0e8, V = -1.2e12, nsteps = 16,
            yield = x -> c + α * ρa * gg * (hh - x[2]), rtol = 1.0e-6)
    end

    # ---- DP WITH A VON MISES CAP (yield-strength envelope) ---- σ_Y(ζ) = min(c + α·ρg·ζ, σ_max).
    # The brittle DP branch makes the BASE too strong (97 MPa) to yield ⇒ only the top (extension) fibre
    # yields. Capping the deep strength at a ductile von-Mises ceiling σ_max lets the COMPRESSION (bottom)
    # fibre yield too ⇒ plastic from BOTH surfaces inward, elastic core between. c (top) and σ_max (bottom)
    # act as two convex, independently-tunable strengths — no Lode convexity cap — so the neutral plane
    # ζn ≈ h·σ_max/(c+σ_max) ≈ 31 km is retained WITH two-sided yielding. min() of two convex surfaces is
    # convex ⇒ stable; the DP-cone/vM-cylinder corner never appears here (σ_Y is a per-point scalar of depth).
    # Load V controls the loaded-edge (side-wall) shear; the cap σ_max controls the base — independent knobs.
    # The loaded face yields when the applied parabolic-shear von Mises √3·1.5·|V|/h exceeds σ_Y(mid); the
    # a-priori limit is |V|_max = σ_Y(mid)·h/(1.5√3) = 1.155 TN/m for σ_Y(mid)=σmax=60. V=−1.0 (the unified
    # project load) sits safely under it: face-shear vM ≈ 52 < 60 ⇒ clean face, while the bulge fibre vM
    # (≈1.6× the face) ≈ 83 > 60 ⇒ BOTH fibres yield (extension top + shortening base), two-sided.
    let c = 40.0e6, σmax = 60.0e6, ζn = 30.0e3
        α = c * (hh - 2ζn) / ((ζn^2 - hh^2 / 2) * ρa * gg)
        run_grav_plastic("dp_capped"; σ₀ = c, H = 7.0e8, V = -1.0e12, nsteps = 16,
            yield = x -> min(c + α * ρa * gg * (hh - x[2]), σmax), rtol = 1.0e-6)
    end

    # ---- STABLE VON MISES (symmetric) — the consolidated working model ----
    # NO end moment (an applied end moment is an EDGE loader → corner hinge/junk). Flexure driven by the
    # shear V, which yields the BULGE (outer rise), not the loaded face. σ_Y is chosen in the edge-elastic
    # window (face_shear, ~1.6·face_shear), face_shear = √3·1.5·|V|/h: keeps the loaded face elastic (no
    # breakdown) while the bulge yields two-sided. Verified: edge 0% yielded, bulge two-sided, stable.
    #   V=1.0 → face 52 MPa, window σ_Y∈(52,83): σ_Y 58–68 (1–25% bulge yield, trench ~0.84 km)
    #   V=1.5 → face 78 MPa, window σ_Y∈(78,125): σ_Y 85–95 (18–29% bulge yield, trench ~1.26 km)  ← default
    run_grav_plastic("vm_sym"; σ₀ = 85.0e6, H = 7.0e8, V = -1.5e12, nsteps = 16, rtol = 1.0e-6)

    # ---- ASYMMETRIC (capped DP/Byerlee) — PAUSED, kept for later (see FINDINGS §9) ----
    # σ_ext(c)=123, σ_comp(cap)=247, R=2 ⇒ neutral plane ~30 km vs the symmetric mid-plate. The version with
    # the geophysical moment M=1e17 hit the loaded-edge corner hinge (moment = edge loader); revisit by
    # driving it with V only (like vm_sym above) and scaling V+σ_Y for amplitude. Left here documented.
    let c = 123.0e6, σmax = 247.0e6, ζn = 30.0e3
        α = c * (hh - 2ζn) / ((ζn^2 - hh^2 / 2) * ρa * gg)
        run_grav_plastic("dp_byerlee"; σ₀ = c, H = 7.0e8, V = -1.0e12, M = 1.0e17, nsteps = 24,
            yield = x -> min(c + α * ρa * gg * (hh - x[2]), σmax), rtol = 1.0e-6)
    end
end
