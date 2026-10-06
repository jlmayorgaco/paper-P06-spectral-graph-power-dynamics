@testset "ExpA bundle serialization round trip" begin
    case = BG.SyntheticSystems.build_synthetic_cases()[3]
    metadata = Dict{String,Any}(
        "source" => "synthetic test",
        "retained_state_names" => ["q$i" for i in 1:5],
        "units" => "SI; s in rad/s",
        "sign_convention" => "self-energy on the left side",
        "operating_point" => "test",
        "model_version" => "test model 1",
        "exact_or_approximate" => "exact",
        "valid_frequency_band" => [0.01, 100.0])
    bundle = BG.ExpAAdapter.bundle_from_realization(case.M, case.L,
        case.nodal_sys, metadata)
    path, io = mktemp()
    close(io)
    try
        BG.ExpAAdapter.save_expA_bundle(path,bundle)
        loaded = BG.ExpAAdapter.load_expA_bundle(path)
        @test loaded.sigma_available
        @test isempty(BG.ExpAAdapter.missing_metadata(loaded))
        s = -0.2 + 2.7im
        @test norm(bundle.sigma(s)-loaded.sigma(s)) < 1e-12
        @test norm(bundle.sigma_derivative(s)-loaded.sigma_derivative(s)) < 1e-12
    finally
        rm(path; force=true)
    end

    unavailable = BG.ExpAAdapter.BNDOperatorBundle(case.M, case.L, nothing, nothing,
        metadata; sigma_available=false)
    path2, io2 = mktemp()
    close(io2)
    try
        BG.ExpAAdapter.save_expA_bundle(path2, unavailable)
        loaded2 = BG.ExpAAdapter.load_expA_bundle(path2)
        @test !loaded2.sigma_available
        @test loaded2.sigma === nothing
    finally
        rm(path2; force=true)
    end
end
