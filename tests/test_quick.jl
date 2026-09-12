# tests/test_quick.jl — fast checks of the constitutive kernel (no solve).   julia --project=. tests/test_quick.jl
# Run from tests/test_quick.py (pytest) or directly.  Exits nonzero on failure.
include(joinpath(@__DIR__, "..", "model", "gpe_plastic.jl"))
using Test

q_inplane(σ) = sqrt((σ[1, 1] - σ[3, 3])^2 + (2σ[1, 3])^2)      # the in-plane Tresca measure used by the return map
comps = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))
unit(i, j) = SymmetricTensor{2, 3}((a, b) -> ((a == i && b == j) || (a == j && b == i)) ? 1.0 : 0.0)

@testset "Tresca return lands on the yield surface" begin
    m = J2Plasticity(70.0e9, 0.25, 150.0e6, 0.0; crit = :tresca); st = PState()
    ϵ = SymmetricTensor{2, 3}((4.0e-3, 0.0, 1.0e-3, 0.0, 0.0, -1.0e-3))       # well beyond yield
    σᵗ = m.Dᵉ ⊡ ϵ
    @test q_inplane(σᵗ) > 150.0e6
    σ, D, st1 = stress_tangent(ϵ, m, st)
    @test isapprox(q_inplane(σ), 150.0e6; rtol = 1.0e-9)                       # on the surface
    @test st1.k > 0                                                            # plastic multiplier accumulated
    @test isapprox(σ[2, 2], σᵗ[2, 2]; rtol = 1.0e-12)                          # out-of-plane component untouched
    # below yield: elastic pass-through, tangent = Dᵉ
    ϵe = 1.0e-4 * ϵ
    σe, De, ste = stress_tangent(ϵe, m, st)
    @test σe ≈ m.Dᵉ ⊡ ϵe && De == m.Dᵉ && ste.k == 0
end

@testset "consistent tangent = finite difference of the return-mapped stress (Tresca and von Mises)" begin
    for crit in (:tresca, :mises)
        m = J2Plasticity(70.0e9, 0.25, 150.0e6, 0.0; crit = crit); st = PState()
        ϵ = SymmetricTensor{2, 3}((4.0e-3, 0.3e-3, 1.0e-3, -0.5e-3, 0.2e-3, -1.0e-3))
        σ, D, _ = stress_tangent(ϵ, m, st)
        @test q_inplane(σ) < q_inplane(m.Dᵉ ⊡ ϵ)                               # plastic regime
        h = 1.0e-9
        for (i, j) in comps
            e = unit(i, j)
            dσ = (stress_tangent(ϵ + h * e, m, st)[1] - stress_tangent(ϵ - h * e, m, st)[1]) / (2h)
            Dad = D ⊡ e
            @test isapprox(dσ, Dad; rtol = 1.0e-5) || (println("crit=$crit comp=($i,$j)\nFD  = $dσ\nAD  = $Dad"); false)
        end
    end
end
