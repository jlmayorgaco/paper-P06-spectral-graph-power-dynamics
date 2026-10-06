@testset "graph dissipative power identity" begin
    cases = BG.SyntheticSystems.build_synthetic_cases()
    rng = MersenneTwister(71)
    errors = Float64[]
    for case in cases
        for f in (0.1, 1.3, 10.0)
            z = randn(rng, ComplexF64, length(case.Lambda))
            result = BG.GraphModalOperator.harmonic_power_identity(case.basis.Phi,
                s -> BG.DynamicSelfEnergy.sigma(case.nodal_sys, s), 2π*f, z)
            push!(errors, result.relative_error)
        end
    end
    @test maximum(errors) < 1e-10
    @test all(isreal(e) for e in errors)
end
