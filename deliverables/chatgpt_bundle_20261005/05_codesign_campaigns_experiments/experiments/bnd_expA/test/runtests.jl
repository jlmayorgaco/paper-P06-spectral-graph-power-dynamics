using LinearAlgebra
using Test

@testset "Experiment A algebra and indexing" begin
    # State order is a bijection over retained, condensed, and algebraic sets.
    retained = [1, 3, 4]
    condensed = [2, 5]
    algebraic = [6, 7]
    allidx = vcat(retained, condensed, algebraic)
    @test length(unique(allidx)) == 7
    @test sort(allidx) == collect(1:7)
    @test size(Matrix{Float64}(I, length(retained), length(retained))) == (3, 3)

    # Synthetic exact Schur identity, evaluated by a solve rather than inverse.
    T = ComplexF64[4+im 1-2im; 3+im 5-im]
    rr, cc = [1], [2]
    S = T[rr,rr] - T[rr,cc] * (T[cc,cc] \ T[cc,rr])
    @test isapprox(det(T), det(T[cc,cc]) * det(S); rtol=1e-13, atol=1e-13)
    lifted = vcat(Matrix{ComplexF64}(I,1,1), -(T[cc,cc] \ T[cc,rr]))
    @test norm(T * lifted - vcat(S, zeros(ComplexF64,1,1))) <= 1e-13

    # The general eigenvalue sensitivity sign for T(λ,z)=λE-A(z).
    E = Matrix{Float64}(I, 2, 2)
    A0 = [0.0 1.0; -2.0 -3.0]
    Az = [0.0 0.0; -1.0 0.0]
    e0 = eigen(A0)
    k = argmin(abs.(e0.values .- (-1.0)))
    λ = e0.values[k]
    V = e0.vectors
    W = inv(V)
    l = conj.(W[k,:])
    r = V[:,k]
    analytic = dot(l, Az*r) / dot(l, E*r)
    h = 1e-6
    λp = eigen(A0 + h*Az).values[argmin(abs.(eigen(A0+h*Az).values .- λ))]
    λm = eigen(A0 - h*Az).values[argmin(abs.(eigen(A0-h*Az).values .- λ))]
    numeric = (λp-λm)/(2h)
    @test isapprox(analytic, numeric; rtol=1e-6, atol=1e-8)
    @test isapprox(analytic, -(dot(l, -Az*r) / dot(l, E*r)); rtol=1e-12)

    # Conjugate-pair and rad/s-to-Hz handling.
    pair = -0.2 + 2pi*0.7im
    @test abs(imag(pair))/(2pi) ≈ 0.7
    @test abs(imag(conj(pair)))/(2pi) ≈ 0.7
    @test abs(0.0)/(2pi) == 0.0
end
