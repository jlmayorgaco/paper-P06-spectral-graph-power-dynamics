using Test, CSV, DataFrames, SHA, TOML, Statistics

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const TABLES=joinpath(ROOT,"reports","experiment_E","tables")
table(name)=CSV.read(joinpath(TABLES,name),DataFrame)

@testset "ExpE frozen analytical evidence" begin
    config=joinpath(ROOT,"experiments","bnd_expE","configs","DESIGN_DOMAIN_FROZEN.toml")
    hash=strip(read(config*".sha256",String))
    @test hash==bytes2hex(sha256(read(config)))
    domain=TOML.parsefile(config)
    @test domain["replaceable_buses"]==collect(30:39)
    @test domain["independent_PLL_gains"]==true
    @test domain["Kp_factor_min"]==domain["Ki_factor_min"]==0.25
    @test domain["Kp_factor_max"]==domain["Ki_factor_max"]==4.0

    generators=table("TABLE_E01_replaceable_generators.csv")
    @test nrow(generators)==10
    @test all(generators.replaceable)
    @test abs(sum(generators.SG_dispatch_initial_MW)-5402.761089978776)<1e-6
    @test abs(only(generators[generators.bus.==39,:SG_dispatch_initial_MW])-271.0422158686306)<1e-6

    loads=table("TABLE_E03_initialized_ZIP_loads.csv")
    @test abs(only(loads[loads.bus.==31,:initialized_load_MW])-0.10780220734605216)<1e-6
    @test abs(only(loads[loads.bus.==39,:initialized_load_MW])-375.04221586863076)<1e-6

    single=table("TABLE_E04_single_bus_spectrum_audit.csv")
    @test sort(single.bus)==[30,33,35,37]
    @test maximum(single.spectral_Hausdorff)<1e-8
    woodbury=table("TABLE_E05_PLL_Woodbury_audit.csv")
    @test median(woodbury.Woodbury_relative_error)<1e-11
    @test quantile(woodbury.Woodbury_relative_error,0.95)<1e-9
    @test maximum(woodbury.A_factor_error)<1e-10
    closure=table("TABLE_E06_collective_closure_audit.csv")
    @test all(closure.closure_dimension.==20)
    @test maximum(closure.logabsdet_error)<1e-10
    @test maximum(closure.phase_error)<1e-10

    neutral=table("TABLE_E07_all_GFL_neutrality.csv")
    @test nrow(neutral)==3
    @test maximum(neutral.rotation_null_residual)<1e-15
    @test all(abs.(neutral.closest_real).<1e-4)
    @test all(abs.(neutral.second_real).<1e-4)

    provisional=table("TABLE_E14_provisional_per_generator_design.csv")
    @test provisional.bus==collect(30:39)
    @test all((0 .<= provisional.rho) .& (provisional.rho .<= 1))
    @test all(provisional.Kp_over_nom .>= 0.25-1e-10)
    @test all(provisional.Kp_over_nom .<= 4+1e-10)
    @test all(provisional.Ki_over_nom .>= 0.25-1e-10)
    @test all(provisional.Ki_over_nom .<= 4+1e-10)
    @test abs(sum(provisional.GFL_MW)-5401.57187809661)<1e-6
    @test abs(sum(provisional.retained_SG_MW)-1.1892118821653361)<1e-6
    balance=table("TABLE_E15_provisional_electrical_balance.csv")
    @test maximum(balance.KCL_abs)<1e-8
    @test !occursin("PowerDynamics",join(read(joinpath(ROOT,"src","bnd_design_e",f),String)
        for f in readdir(joinpath(ROOT,"src","bnd_design_e")) if endswith(f,".jl")))
end
