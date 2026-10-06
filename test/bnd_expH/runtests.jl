using Test, CSV, DataFrames, TOML, SHA

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPORT = joinpath(ROOT, "reports", "experiment_H")
const TABLES = joinpath(REPORT, "tables")

@testset "ExpH frozen analytical and post-freeze evidence" begin
    integrity = CSV.read(joinpath(TABLES, "TABLE_H00_model_integrity.csv"), DataFrame)
    @test all(.!occursin.("FAIL", integrity.status))
    @test any(occursin.("PASS_ARCHIVED_FULL_SPECTRUM", integrity.status))

    quotient = CSV.read(joinpath(TABLES, "TABLE_H01_gauge_quotient.csv"), DataFrame)
    @test nrow(quotient) == 6
    @test all(quotient.classification .== "SIMPLE_ZERO_AFTER_QUOTIENT")
    @test all(quotient.quotient_nullity .== 1)
    @test all(quotient.quotient_algebraic_multiplicity .== 1)
    @test all(quotient.full_algebraic_multiplicity .== 2)

    scaling = CSV.read(joinpath(TABLES, "TABLE_H03_scaling_exponent.csv"), DataFrame)
    @test all(scaling.full_range_within_5pct)
    @test all(abs.(scaling.full_range_exponent .- 1) .<= 0.05)

    derivative = CSV.read(joinpath(TABLES, "TABLE_H05_derivative_validation.csv"), DataFrame)
    @test all(derivative[derivative.metric .!= "near_zero_cutoff", :status][1:2] .== "PASS")
    rank = CSV.read(joinpath(TABLES, "TABLE_H05_rank.csv"), DataFrame)
    @test all(rank.numerical_rank .== 1)
    active = CSV.read(joinpath(TABLES, "TABLE_H06_active_basis.csv"), DataFrame)
    @test last(active.cumulative_energy) >= 0.999

    final_path = joinpath(REPORT, "Z_H_FINAL.toml")
    candidate = TOML.parsefile(final_path)
    sha = bytes2hex(sha256(read(final_path)))
    @test strip(read(final_path * ".sha256", String)) == sha
    @test candidate["candidate_frozen"]
    @test candidate["design_used_PowerDynamics"] == false
    @test candidate["beta_required"] == 1e-4
    @test isapprox(candidate["retained_SG_MW"], 11.2; atol=1e-10)
    @test isapprox(candidate["spectral_abscissa_s_inv"], -0.06479596417248183; atol=1e-12)

    frontier = CSV.read(joinpath(TABLES, "TABLE_H14_robustness_frontier.csv"), DataFrame)
    robust = only(frontier[frontier.beta_required .== 1e-4, :])
    @test robust.robust_pass
    @test robust.spectral_abscissa <= -0.05
    @test robust.small_gain_margin > 0
    branches = CSV.read(joinpath(TABLES, "TABLE_H16_branch_summary.csv"), DataFrame)
    @test any(branches.robust_feasible)
    kkt = CSV.read(joinpath(TABLES, "TABLE_H15_KKT_iterations.csv"), DataFrame)
    @test only(kkt.iterations) == 0
    @test occursin("NOT_SOLVED", only(kkt.status)) || occursin("NOT_CERTIFIED", only(kkt.status))

    design_source = read(joinpath(ROOT, "src", "bnd_design_h", "BNDDesignH.jl"), String)
    runner_source = read(joinpath(ROOT, "experiments", "bnd_expH", "run_experiment_H.jl"), String)
    @test !occursin(r"(?im)^\s*(using|import)\s+PowerDynamics", design_source)
    @test !occursin(r"(?im)^\s*(using|import)\s+PowerDynamics", runner_source)

    pd = CSV.read(joinpath(TABLES, "TABLE_H20_powerdynamics_validation.csv"), DataFrame)
    @test only(pd.candidate_sha256) == sha
    @test only(pd.status) == "EVALUATED"
    @test only(pd.candidate_fixed_point)
    @test only(pd.candidate_equilibrium_residual_inf) <= 1e-9
    @test only(pd.alpha_PD_s_inv) > -0.05
    @test !only(pd.meets_strict_margin)

    tds = CSV.read(joinpath(TABLES, "TABLE_H21_tds_validation.csv"), DataFrame)
    @test nrow(tds) == 2
    @test all(tds.status .== "EVALUATED")
    @test all(tds.finite)
    @test maximum(tds.relative_frequency_scaling_error) <= 0.10
    @test maximum(tds.relative_voltage_scaling_error) <= 0.10
end
