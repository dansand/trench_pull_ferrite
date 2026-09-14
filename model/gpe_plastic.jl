# gpe_plastic.jl — the finite-element solver behind every model in the paper: a total-Lagrangian, plane-strain,
# elasto-plastic plate on a Winkler (Δρg) foundation, loaded by a parabolic end shear at the trench face.
#
# Small-strain J2 / in-plane Tresca plasticity in the MATERIAL frame: the return map acts on the Green–Lagrange
# strain E and the 2nd-Piola stress S (both rotation-invariant), so the large flexural rotation is carried exactly by
# the total-Lagrangian kinematics and — since |E| is small — S ≈ σ_Cauchy (co-rotational J2 to O(strain)).
# Perfect plasticity (H = 0) for the uniform Tresca models; linear hardening H with a depth-dependent von Mises
# yield strength for the DD-VM models.  Driven by paper_models.jl (production suites) and idealized_beam.jl (the
# SI benchmarks); never run on its own.
#
# 2026-09-14: stripped to the code paths the shipped models exercise (follower foundation, displacement-controlled
# loading, end moment, symmetry-plane loading, Lode-angle strength, the "bending" return and the pre-plasticity
# elastic track gpe_mvm.jl removed; regeneration verified byte-identical afterwards).

using Ferrite, Tensors, LinearAlgebra, Printf

const ẑ = Vec{2}((0.0, 1.0))

# ---- J2 plane-strain return map ----
struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}}
    G::T; K::T; σ₀::T; H::T; crit::Symbol; Dᵉ::S
end
function J2Plasticity(E, ν, σ₀, H; crit = :mises)
    δ(i, j) = i == j ? 1.0 : 0.0
    G = E / 2(1 + ν); K = E / 3(1 - 2ν)
    f(i, j, k, l) = 2.0G * (0.5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) + ν / (1.0 - 2.0ν) * δ(i, j) * δ(k, l))
    return J2Plasticity(G, K, σ₀, H, crit, SymmetricTensor{4, 3}(f))
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
    # Δϵᵖ = Dᵉ⁻¹:(σᵗ−σ) reduces EXACTLY to (σᵗ−σ)/2G: the in-plane return keeps (σxx+σzz)/2 and σyy, so tr(σᵗ−σ) = 0
    # and the deviatoric compliance is 1/2G (tests/test_quick.jl checks the two agree to machine precision).
    return σ, D, PState(st.ϵᵖ + (σᵗ - σ) / (2 * m.G), σ, st.k + μ)
end

# trial → yield check → radial return → consistent tangent (small strain, 3D).  σ0_local carries the
# depth-dependent von Mises strength of the DD-VM models (constant = m.σ₀ otherwise).
function stress_tangent(ϵ::SymmetricTensor{2, 3}, m::J2Plasticity, st::PState, σ0_local = m.σ₀)
    m.crit === :tresca && return tresca_return(ϵ, m, st, σ0_local)
    G, H = m.G, m.H
    σᵗ = m.Dᵉ ⊡ (ϵ - st.ϵᵖ); sᵗ = dev(σᵗ); σᵗₑ = sqrt(1.5 * sᵗ ⊡ sᵗ)
    σy = σ0_local
    φ = σᵗₑ - (σy + H * st.k)
    if φ < 0.0
        return σᵗ, m.Dᵉ, PState(st.ϵᵖ, σᵗ, st.k)
    end
    h = H + 3G; μ = φ / h                    # φ ≥ 0 here (the φ < 0 branch returned above) ⇒ μ ≥ 0 without a clamp
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
# bfn = body force; s0fn = ISOTROPIC dead prestress (added to P, no geometric tangent).  Both are zero in every
# production model (massless platform); the arguments are kept so the assembly is unchanged.
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

function assemble_plastic!(K, g, dh, cv, fv, m, u, states, states_old, springs, tractfn, tractset, bfn, s0fn, yieldfn)
    asm = start_assemble(K, g); nc = getnbasefunctions(cv)
    ke = zeros(nc, nc); ge = zeros(nc)
    for cell in CellIterator(dh)
        c = cellid(cell); reinit!(cv, cell)
        cell_plastic!(ke, ge, cv, m, u[celldofs(cell)], view(states, :, c), view(states_old, :, c), bfn, s0fn, yieldfn, getcoordinates(cell))
        assemble!(asm, celldofs(cell), ke, ge)
    end
    nf = getnbasefunctions(fv); kef = zeros(nf, nf); gef = zeros(nf)
    for (sset, kf, fref) in springs                         # foundation springs: vertical Winkler dead load kf·uz − fref
        for fc in FacetIterator(dh, sset)
            reinit!(fv, fc); fill!(kef, 0); fill!(gef, 0); ue = u[celldofs(fc)]
            for qp in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, qp); uz = function_value(fv, qp, ue) ⋅ ẑ
                for i in 1:nf
                    di = shape_value(fv, qp, i) ⋅ ẑ
                    gef[i] += (kf * uz - fref) * di * dΓ
                    for j in 1:nf; kef[i, j] += kf * di * (shape_value(fv, qp, j) ⋅ ẑ) * dΓ; end
                end
            end
            assemble!(asm, celldofs(fc), kef, gef)
        end
    end
    for fc in FacetIterator(dh, tractset)                   # end shear traction (+ the confining/membrane σxx)
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
function solve_plastic(; L, h, nx, nz, E, ν, σ₀, H, springs, tract_z,
                         bodyforce = x -> Vec{2}((0.0, 0.0)), prestress = x -> ZERO_S0₂, confine_x = x -> 0.0,
                         yield = nothing, crit = :mises, order = 2, nsteps = 10, rtol = 1.0e-8,
                         load_face = "right", clamp_face = "left")
    m = J2Plasticity(E, ν, σ₀, H; crit = crit)
    yieldfn = yield === nothing ? (x -> σ₀) : yield        # depth-dependent von Mises: σ_Y(x); default = constant σ₀
    grid = generate_grid(Quadrilateral, (nx, nz), Vec(0.0, 0.0), Vec(L, h))
    ip = Lagrange{RefQuadrilateral, order}()^2; qo = 2 * order
    cv = CellValues(QuadratureRule{RefQuadrilateral}(qo), ip)
    fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(qo), ip)
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    ch = ConstraintHandler(dh); add!(ch, Dirichlet(:u, getfacetset(grid, clamp_face), (x, t) -> [0.0, 0.0], [1, 2]))
    close!(ch)
    K = allocate_matrix(dh); u = zeros(ndofs(dh)); g = zeros(ndofs(dh))
    nqp = getnquadpoints(cv); ncells = getncells(grid)
    states = [PState() for _ in 1:nqp, _ in 1:ncells]; states_old = [PState() for _ in 1:nqp, _ in 1:ncells]
    springsets = [(getfacetset(grid, nm), kf, fr) for (nm, kf, fr) in springs]
    tractset = getfacetset(grid, load_face); hmid = h / 2
    for step in 1:nsteps
        f = step / nsteps
        bfn = x -> f * bodyforce(x); s0fn = x -> f * prestress(x)
        shfn = z -> 1.5 * tract_z * (1 - ((z - hmid) / hmid)^2)                 # St-Venant parabola, 0 at both surfaces
        tractfn = x -> Vec{2}((f * confine_x(x), f * shfn(x[2])))               # trench face: (confine, shear)
        springs_f = [(set, kf, f * fr) for (set, kf, fr) in springsets]
        update!(ch, f); apply!(u, ch); ref = 1.0
        for it in 1:40
            assemble_plastic!(K, g, dh, cv, fv, m, u, states, states_old, springs_f, tractfn, tractset, bfn, s0fn, yieldfn)
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

# Top-surface topography: u_z(x) of the nodes on the reference surface y = y_surf, sorted by x.  Native z-down mesh:
# surface at y = 0 (pass y_surf = 0.0) and w is then positive-DOWN.  Uses evaluate_at_grid_nodes (the global dof index
# is NOT 2i per grid node).
function topography(res, h; y_surf = h)
    uv = Ferrite.evaluate_at_grid_nodes(res.dh, res.u, :u)
    nodes = [res.grid.nodes[i].x for i in 1:getnnodes(res.grid)]
    top = [i for i in eachindex(nodes) if abs(nodes[i][2] - y_surf) < 1.0e-6]
    x = [nodes[i][1] for i in top]; w = [uv[i][2] for i in top]
    o = sortperm(x)
    return x[o], w[o]
end

# export the 2nd-PK stress (+ the isotropic prestress, zero here), von Mises, accumulated plastic strain and the yield
# flag from the committed states, plus the consistent Cauchy shear and the FE spatial gradients the Python side needs.
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
