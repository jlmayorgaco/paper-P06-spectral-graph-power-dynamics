@testset "exact Schur self-energy and determinant factorization" begin
    cases = BG.SyntheticSystems.build_synthetic_cases()
    case = cases[3]
    sys = BG.DynamicSelfEnergy.with_offdiagonal_scale(case.modal_sys, 0.4)
    s = -0.15 + 1im * 2π * 1.2
    S = BG.DynamicSelfEnergy.sigma(sys, s)
    T = s^2 * I + s * S + Diagonal(case.Lambda)
    for k in 2:length(case.Lambda)
        gamma = BG.IntermodalSelfEnergy.schur_self_energy(T, k)
        @test gamma.residual < 1e-13
        @test BG.IntermodalSelfEnergy.determinant_factorization_error(T, k) < 1e-10
        @test isapprox(gamma.effective, T[k,k] + gamma.gamma; rtol=1e-13, atol=1e-13)
    end
    s0 = BG.SyntheticSystems.build_synthetic_cases()[1]
    T0 = s^2 * I + s * BG.DynamicSelfEnergy.sigma(s0.modal_sys, s) + Diagonal(s0.Lambda)
    @test maximum(abs(BG.IntermodalSelfEnergy.schur_self_energy(T0,k).gamma) for k in 2:5) < 1e-12
end
