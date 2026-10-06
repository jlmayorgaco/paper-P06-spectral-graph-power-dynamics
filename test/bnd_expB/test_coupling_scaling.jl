@testset "quadratic coupling scaling and pairwise reconstruction" begin
    case = BG.SyntheticSystems.build_synthetic_cases()[4]
    k = 2
    s = -0.12 + 1im * 2π * 0.8
    epsilons = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
    gammas = Float64[]
    for e in epsilons
        sys = BG.DynamicSelfEnergy.with_offdiagonal_scale(case.modal_sys, e)
        S = BG.DynamicSelfEnergy.sigma(sys, s)
        T = s^2 * I + s * S + Diagonal(case.Lambda)
        push!(gammas, abs(BG.IntermodalSelfEnergy.schur_self_energy(T,k).gamma))
    end
    slope = (log(last(gammas)) - log(first(gammas))) /
            (log(last(epsilons)) - log(first(epsilons)))
    @test 1.9 <= slope <= 2.1

    eps0 = first(epsilons)
    baseS = BG.DynamicSelfEnergy.sigma(case.modal_sys, s)
    approx = eps0^2 * BG.IntermodalSelfEnergy.pairwise_second_order(baseS, case.Lambda, s, k)
    sys = BG.DynamicSelfEnergy.with_offdiagonal_scale(case.modal_sys, eps0)
    T = s^2 * I + s * BG.DynamicSelfEnergy.sigma(sys, s) + Diagonal(case.Lambda)
    exact = BG.IntermodalSelfEnergy.schur_self_energy(T,k).gamma
    @test abs(exact - approx) / max(abs(exact), eps(Float64)) < 0.01
end
