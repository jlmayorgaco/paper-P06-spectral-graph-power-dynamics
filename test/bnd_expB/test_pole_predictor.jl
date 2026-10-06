@testset "pole predictor, tracking, and augmented-system ground truth" begin
    case = BG.SyntheticSystems.build_synthetic_cases()[4]
    k = 2
    s0 = BG.PolePredictor.uncoupled_modal_pole(case.modal_sys, case.Lambda, k;
        mass=Matrix{Float64}(I,5,5))
    eps0 = 1e-4
    sys = BG.DynamicSelfEnergy.with_offdiagonal_scale(case.modal_sys, eps0)
    prediction = BG.PolePredictor.predict_pole_shift(sys, case.Lambda, s0, k)
    A = BG.DynamicSelfEnergy.augmented_matrix(Matrix{Float64}(I,5,5),
        Diagonal(case.Lambda), sys)
    exact_pair = BG.Reporting.positive_imag_eigenpair(A, s0)
    @test abs(exact_pair.value - prediction.predicted) / abs(exact_pair.value) < 0.01
    @test abs(prediction.shift) > 0
    @test abs(prediction.gamma) > 0

    # Real augmented matrices must retain conjugate pairs.
    vals = eigvals(A)
    for v in filter(v -> imag(v) > 1e-7, vals)
        @test minimum(abs.(vals .- conj(v))) < 1e-9
    end

    # A full-system eigenvalue is a zero of the exact graph-transformed NEP.
    s2 = BG.SyntheticSystems.build_synthetic_cases()[3]
    A2 = BG.DynamicSelfEnergy.augmented_matrix(s2.M, s2.L, s2.nodal_sys)
    eig2 = eigvals(A2)
    target_f = s2.modal_frequencies_hz[2]
    roots = filter(v -> imag(v) > 1e-8, eig2)
    root = roots[argmin(abs.(imag.(roots) .- 2π*target_f))]
    T = BG.GraphModalOperator.graph_operator(s2.basis,
        s -> BG.DynamicSelfEnergy.sigma(s2.nodal_sys,s), root)
    @test minimum(svdvals(T)) / max(1.0, opnorm(T,2)) < 1e-8

    q = randn(MersenneTwister(9), ComplexF64, length(s2.Lambda))
    x = (root * I - s2.nodal_sys.Ac) \ (root * s2.nodal_sys.B * q)
    direct = root^2 * s2.M * q + root * s2.nodal_sys.D0 * q +
             s2.L * q + s2.nodal_sys.C * x
    reduced = (root^2 * s2.M + root * BG.DynamicSelfEnergy.sigma(s2.nodal_sys,root) + s2.L) * q
    @test norm(direct-reduced) / max(1.0,norm(direct)) < 1e-10

    # Analytic derivative agrees with a centered complex finite-difference check.
    s = -0.3 + 4.1im
    h = 1e-6
    analytic = BG.DynamicSelfEnergy.sigma_derivative(s2.modal_sys,s)
    finite = (BG.DynamicSelfEnergy.sigma(s2.modal_sys,s+h) -
              BG.DynamicSelfEnergy.sigma(s2.modal_sys,s-h)) / (2h)
    @test norm(analytic-finite) / max(1.0,norm(analytic)) < 1e-8
end
