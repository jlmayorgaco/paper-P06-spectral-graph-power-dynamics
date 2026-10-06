@testset "generalized graph basis" begin
    R = [1.2 0.1 0.03 0.0; 0.0 0.9 0.06 0.02; 0.0 0.0 1.1 0.04; 0.0 0.0 0.0 1.3]
    M = R' * R
    L = R' * Diagonal([0.0, 2.0, 2.0, 9.0]) * R
    basis = BG.GraphBasis.generalized_graph_basis(M, L; tol_deg=1e-9)
    @test basis.M_orthogonality_error < 1e-12
    @test basis.L_diagonalization_error < 1e-10
    @test length(basis.zero_modes) == 1
    @test BG.GraphBasis.oscillatory_indices(basis) == [2, 3, 4]
    @test any(length(g) == 2 for g in basis.clusters)

    cluster = [2, 3]
    A = [1.0 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
    U = Matrix{Float64}(I, 4, 4)
    theta = 0.41
    U[cluster, cluster] = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    @test isapprox(BG.GraphBasis.block_frobenius(A, cluster),
                   BG.GraphBasis.block_frobenius(U' * A * U, cluster); rtol=1e-14)
    @test_throws PosDefException BG.GraphBasis.generalized_graph_basis(zeros(2,2), Matrix{Float64}(I,2,2))

    cases = BG.SyntheticSystems.build_synthetic_cases()
    @test cases[5].modal_frequencies_hz == [0.0, 3.0, 3.005, 6.0, 9.0]
    case = cases[3]
    s = -0.2 + 3.7im
    Tnode = s^2 * case.M + s * BG.DynamicSelfEnergy.sigma(case.nodal_sys, s) + case.L
    Tgraph = case.basis.Phi' * Tnode * case.basis.Phi
    Sgraph = BG.GraphModalOperator.graph_sigma(case.basis,
        x -> BG.DynamicSelfEnergy.sigma(case.nodal_sys, x), s)
    Texpected = s^2 * I + s * Sgraph + Diagonal(case.Lambda)
    @test norm(Tgraph - Texpected) / max(1.0, norm(Tgraph)) < 1e-11
end
