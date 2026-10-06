using Test, CSV, DataFrames, LinearAlgebra, SHA, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPORT = joinpath(ROOT, "reports", "experiment_G")
const TABLES = joinpath(REPORT, "tables")

include(joinpath(ROOT, "src", "bnd_design_g", "BNDDesignG.jl"))
using .BNDDesignG

@testset "ExpG frozen candidate and analytical utilities" begin
    final_path = joinpath(REPORT, "Z_G_FINAL.toml")
    sidecar = strip(read(final_path * ".sha256", String))
    @test sidecar == bytes2hex(sha256(read(final_path)))
    candidate = TOML.parsefile(final_path)
    @test candidate["status"] == "FROZEN"
    @test candidate["candidate_frozen"] == true
    @test candidate["design_used_PowerDynamics"] == false
    @test isfile(joinpath(REPORT, "Z_G_CANDIDATE_PREBLIND.toml"))
    @test isfile(joinpath(REPORT, "Z_G_CANDIDATE_PREBLIND.toml.sha256"))

    generators = CSV.read(joinpath(TABLES, "TABLE_G13_final_per_generator_design.csv"), DataFrame)
    @test generators.bus == collect(30:39)
    @test all((0 .<= generators.rho_converted) .& (generators.rho_converted .<= 1))
    @test all((0 .<= generators.epsilon_retained) .& (generators.epsilon_retained .<= 1))
    @test abs(sum(generators.P_i_MW) - 5402.761089978776) < 1e-6
    @test abs(sum(generators.retained_SG_MW) - candidate["retained_SG_MW"]) < 1e-8
    @test abs(sum(generators.converted_GFL_MW) - candidate["GFL_MW"]) < 1e-8

    A = [-1.0 0.0; 0.0 -2.0]
    cert = BNDDesignG.robustness_certificate(A, 0.1, 0.1)
    @test cert.small_gain_pass
    @test cert.nominal_abscissa == -1.0
    @test cert.beta_star > 0.1

    lp = BNDDesignG.enumerate_surrogate([1.0, 2.0, 3.0], [1.0, 2.0, 1.0],
        [1.0, 0.5, 3.0], 1.5, 1.0)
    @test lp.feasible
    @test all(-1e-10 .<= lp.eps .<= 1 + 1e-10)
    @test dot([1.0, 2.0, 1.0], lp.eps) >= 1.5 - 1e-9
    @test dot([1.0, 0.5, 3.0], lp.eps) >= 1.0 - 1e-9

    J = [1.0 0.0 0.0; 0.0 1.0 1.0]
    correction = BNDDesignG.gain_correction(J, [0.4, 0.4], [0.5, 0.5, 0.5],
        zeros(3), ones(3))
    @test norm(correction.residual_after) < 1e-10
    @test all(0 .<= [0.5, 0.5, 0.5] + correction.delta .<= 1)

    direct = CSV.read(joinpath(TABLES, "TABLE_G11_direct_modal_validation.csv"), DataFrame)
    @test "status" ∉ names(direct)
    @test maximum(abs.(direct.rocof_direct .- direct.rocof_modal)) /
        max(maximum(abs.(direct.rocof_direct)), 1e-12) < 1e-8
    @test maximum(abs.(direct.frequency_direct .- direct.frequency_modal)) /
        max(maximum(abs.(direct.frequency_direct)), 1e-12) < 1e-6
    capacity = CSV.read(joinpath(TABLES, "TABLE_G12_transient_capacity.csv"), DataFrame)
    @test only(capacity.status) == "EVALUATED"
    @test isfinite(only(capacity.deltaP_max_total_MW))

    pd = CSV.read(joinpath(TABLES, "TABLE_G16_powerdynamics_spectral_validation.csv"), DataFrame)
    pdc = only(pd[pd.case .== "ExpG_candidate", :])
    @test pdc.status == "EVALUATED"
    @test pdc.fixed_point
    @test pdc.equilibrium_residual_inf < 1e-8
    @test pdc.meets_strict_margin
    @test abs(pdc.delta_alpha_PD_minus_ExpG_s_inv) < 1e-3
    tds = CSV.read(joinpath(TABLES, "TABLE_G17_powerdynamics_nonlinear_tds.csv"), DataFrame)
    @test nrow(tds) == 2
    @test all(tds.status .== "EVALUATED")
    @test all(tds.finite)
    @test all(isfinite, tds.max_abs_rocof_Hz_s)
    @test maximum(tds.relative_frequency_scaling_error) < 0.01
    @test maximum(tds.relative_voltage_scaling_error) < 0.01
    ledger = CSV.read(joinpath(TABLES, "TABLE_G00_acceptance_gate_ledger.csv"), DataFrame)
    @test only(ledger[startswith.(ledger.gate, "I independent PowerDynamics"), :status]) == "PASS"
    @test only(ledger[startswith.(ledger.gate, "J nonlinear TDS"), :status]) == "PASS"

    design_sources = join(read(joinpath(ROOT, "src", "bnd_design_g", f), String)
        for f in readdir(joinpath(ROOT, "src", "bnd_design_g")) if endswith(f, ".jl"))
    @test !occursin("using PowerDynamics", design_sources)
    @test !occursin("import PowerDynamics", design_sources)
end
