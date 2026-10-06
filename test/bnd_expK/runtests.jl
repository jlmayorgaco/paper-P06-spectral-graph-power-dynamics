using Test, LinearAlgebra, CSV, DataFrames, TOML, SHA
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK

@testset "ExpK fixed-support analytical integrity" begin
    ctx=design_context(ROOT)
    @test ctx.buses==collect(30:39)
    @test isapprox(sum(ctx.power),5402.761089978776;atol=1e-8)
    @test ctx.net.no_load_kcl_inf<1e-9
    for b in (31,39)
        row=only(eachrow(ctx.net.load_audit[ctx.net.load_audit.bus .== b,:]))
        @test abs(row.initialized_load_MW-row.CSV_setpoint_load_MW)>1e-3
    end
    @test length(unique(support_mask(support_epsilon(ctx,m,0.5)) for m in 0:1023))==1024
    gfl=evaluate(ctx,zeros(10))
    sg=evaluate(ctx,ones(10))
    @test gfl.model.n_dynamic==90
    @test sg.model.n_dynamic==114
    @test length(gfl.lambda)==89
    @test length(sg.lambda)==113
    @test abs(gfl.lambda[argmin(abs.(gfl.lambda))])<1e-7
    @test abs(sg.alpha-(-0.09806540933624713))<1e-8
    @test gfl.gauge_residual<1e-12
    @test sg.gauge_residual<1e-12
    for (i,b) in enumerate(ctx.buses)
        e=support_epsilon(ctx,1<<(i-1),1e-8)
        v=evaluate(ctx,e)
        @test v.model.n_dynamic==(b==39 ? 96 : 102)
        @test v.gauge_residual<1e-12
        @test abs(v.alpha-maximum(real.(v.lambda)))<1e-12
    end
    v=evaluate(ctx,fill(0.5,10);vectors=true)
    j=first(v.active)
    g=fixed_support_gradient(ctx,v,j)
    @test all(isfinite,real.(g.epsilon))
    @test all(isfinite,real.(g.kp))
    @test all(isfinite,real.(g.ki))
    for group in (:epsilon,:kp,:ki), i in (1,9)
        d=derivative_check(ctx,v,j,group,i;h=group==:epsilon ? 1e-4 : 1e-3)
        @test d.absolute_error<1e-5
    end
    x=[0.0,0.0];lo=[-1.0,-2.0];hi=[2.0,3.0]
    s=gain_box_support([4.0,-5.0],x,lo,hi)
    @test s.delta==[2.0,-2.0]
    @test s.value==18.0
    c=minimum_gain_correction([1.0,1.0],-1.0,x,lo,hi,[1.0,1.0])
    @test c.feasible && dot([1.0,1.0],c.delta)<=-1+1e-10
end

@testset "ExpK generated ledgers" begin
    tab=joinpath(ROOT,"reports","experiment_K","tables")
    if isfile(joinpath(tab,"TABLE_K01_support_catalog.csv"))
        d=CSV.read(joinpath(tab,"TABLE_K01_support_catalog.csv"),DataFrame)
        @test nrow(d)==1024
        @test sort(Int.(d.support_mask))==collect(0:1023)
        @test all(d.gauge_residual .< 1e-8)
        @test !any(d.infeasibility_certified)
    end
    if isfile(joinpath(ROOT,"reports","experiment_K","Z_K_NOMINAL_FINAL.toml"))
        path=joinpath(ROOT,"reports","experiment_K","Z_K_NOMINAL_FINAL.toml")
        @test strip(read(path*".sha256",String))==bytes2hex(sha256(read(path)))
        frozen=TOML.parsefile(path)
        @test frozen["global_optimality_certified"]==false
        @test frozen["design_used_PowerDynamics"]==false
        ctx=design_context(ROOT)
        e=CSV.read(joinpath(ROOT,"reports","experiment_E","tables",
            "TABLE_E14_provisional_per_generator_design.csv"),DataFrame)
        g=CSV.read(joinpath(ROOT,"reports","experiment_G","tables",
            "TABLE_G13_final_per_generator_design.csv"),DataFrame)
        sort!(e,:bus);sort!(g,:bus)
        ee=evaluate(ctx,1 .- Float64.(e.rho),Float64.(e.Kp),Float64.(e.Ki))
        gg=evaluate(ctx,Float64.(g.epsilon_retained),Float64.(g.Kp),Float64.(g.Ki))
        @test abs(ee.alpha-(-0.050000000460885985))<1e-8
        @test abs(gg.alpha-(-0.13793268780355464))<1e-8
        @test ee.retained_mw<frozen["retained_SG_MW"]
    end
end
