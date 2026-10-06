using CSV, DataFrames, LinearAlgebra, Random, Statistics, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F4");mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
const N=ExpP.PDExactDesignN
function main()
    ctx=N.design_context(ROOT);rng=MersenneTwister(20261003)
    rows=NamedTuple[]
    for j in 1:200
        count=rand(rng,1:10);perm=randperm(rng,10);active=perm[1:count]
        epsv=zeros(10);epsv[active].=0.01 .+0.99.*rand(rng,count)
        epsv[active[argmax(epsv[active])]]=max(epsv[active[argmax(epsv[active])]],0.5)
        rho=1 .-epsv
        kp=ctx.kpmin .+rand(rng,10).*(ctx.kpmax.-ctx.kpmin)
        ki=ctx.kimin .+rand(rng,10).*(ctx.kimax.-ctx.kimin)
        m=N.descriptor(ctx,rho,kp,ki)
        red=LinearSecurity.reduced_step_model(ctx,m,rho,16;gauge_vector=N.gauge_vector)
        h0=real.(red.C_bus*(-(red.A\red.B)).+red.D_bus)
        signed=mean(h0)
        keff=signed<0 ? -1/signed : NaN
        push!(rows,(;design=j,support_buses=join((29 .+ findall(epsv.>0)),";"),
            epsilon=join(round.(epsv,digits=9),";"),Kp=join(round.(kp,digits=5),";"),
            Ki=join(round.(ki,digits=5),";"),signed_H0_mean_Hz_per_MW=signed,
            H0_bus_spread_Hz_per_MW=maximum(h0)-minimum(h0),Keff_MW_per_Hz=keff,
            conditioning_A=cond(red.A)))
        println("F4_POINT=",j," support_size=",count," Keff=",keff)
    end
    tbl=DataFrame(rows);CSV.write(joinpath(OUT,"TABLE_Q2_F4_steady_frequency_samples.csv"),tbl)
    good=filter(r->isfinite(r.Keff_MW_per_Hz),eachrow(tbl))
    train=collect(good[1:100]);test=collect(good[101:end])
    X(r)=[1.0; parse.(Float64,split(r.epsilon,";"))]
    A=hcat((X(r) for r in train)... )'
    y=Float64[r.Keff_MW_per_Hz for r in train]
    coef=A\y
    errs=Float64[];rels=Float64[]
    for r in test
        pred=dot(coef,X(r));e=pred-r.Keff_MW_per_Hz
        push!(errs,abs(e));push!(rels,abs(e)/max(abs(r.Keff_MW_per_Hz),1e-12))
    end
    maxrel=maximum(rels);maxabs=maximum(errs)
    status=maxrel<=1e-8 ? "EXACT_AFFINE_STIFFNESS" : maxrel<=1e-3 ? "APPROXIMATE_AFFINE_STIFFNESS" : "NONAFFINE"
    CSV.write(joinpath(OUT,"TABLE_Q2_F4_affine_holdout.csv"),DataFrame(
        holdout_design=[r.design for r in test],Keff_actual=[r.Keff_MW_per_Hz for r in test],
        Keff_predicted=[dot(coef,X(r)) for r in test],absolute_error=errs,relative_error=rels))
    result=Dict("stage"=>"F4","withheld_designs"=>nrow(tbl),"finite_holdout_designs"=>length(test),
        "STEADY_FREQ_STRUCTURE"=>status,"K_load_MW_per_Hz"=>coef[1],
        "c_i_MW_per_Hz"=>coef[2:end],"affine_holdout_max_abs_error_MW_per_Hz"=>maxabs,
        "affine_holdout_max_relative_error"=>maxrel,"fit_training_designs"=>length(train),
        "finite_steady_frequency_outputs"=>length(good),"model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a")
    open(joinpath(OUT,"F4_RESULTS.toml"),"w") do io;TOML.print(io,result);end
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
    # F4 — Steady-frequency structure

    Tested 200 deterministic mixed/support designs. For each, `H0=-C A⁻¹B+D` was computed from the exact reduced model at all generator-bus frequency outputs. The scalar stiffness fit `Keff=-1/mean(H0)` uses 100 training cases and 100 holdout cases; input gains were sampled over the existing frozen bounds. This is a diagnostic affine fit, not a global certificate.

    - Classification: **$(status)**.
    - Holdout max relative error: $(maxrel); max absolute error $(maxabs) MW/Hz.
    - Fitted Kload: $(coef[1]) MW/Hz.
    - Fitted c_i: $(join(coef[2:end],", ")).
    """)
    println("F4_STATUS=",status," max_relative_holdout_error=",maxrel)
end
main()
