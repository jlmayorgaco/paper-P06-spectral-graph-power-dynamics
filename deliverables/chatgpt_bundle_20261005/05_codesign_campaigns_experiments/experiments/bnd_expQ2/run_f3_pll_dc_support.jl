using CSV, DataFrames, LinearAlgebra, Random, Statistics, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F3");mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(ROOT,"src","bnd_expQ2","GridFrequency.jl"))
include(joinpath(ROOT,"src","bnd_design_g","Robustness.jl"))
const N=ExpP.PDExactDesignN

function model_dc(ctx,rho,kp,ki)
    m=N.descriptor(ctx,rho,kp,ki)
    red=LinearSecurity.reduced_step_model(ctx,m,rho,16;gauge_vector=N.gauge_vector)
    h0=real.(red.C_bus*(-(red.A\red.B)).+red.D_bus)
    alpha=maximum(real.(eigvals(red.A)))
    (;h0,alpha,m,red)
end
function beta_at_frequency(A,omega,sigma=0.05)
    R=(im*omega*I-(A+sigma*I))\I
    F=svd(R)
    (;beta=1/F.S[1],sigma_peak=F.S[1])
end

function main()
    ctx=N.design_context(ROOT);rng=MersenneTwister(20261002)
    rows=NamedTuple[];maxkp=0.0;maxki=0.0;maxalpha=0.0
    for j in 1:30
        rho=0.15 .+ 0.8 .* rand(rng,10)
        kp=ctx.kpmin .+ rand(rng,10).*(ctx.kpmax.-ctx.kpmin)
        ki=ctx.kimin .+ rand(rng,10).*(ctx.kimax.-ctx.kimin)
        busidx=mod1(j,10);gain=iseven(j) ? "Kp" : "Ki"
        base=model_dc(ctx,rho,kp,ki)
        Δ=gain=="Kp" ? 1e-4*(ctx.kpmax[busidx]-ctx.kpmin[busidx]) :
                       1e-4*(ctx.kimax[busidx]-ctx.kimin[busidx])
        kpP=copy(kp);kpM=copy(kp);kiP=copy(ki);kiM=copy(ki)
        gain=="Kp" ? (kpP[busidx]+=Δ;kpM[busidx]-=Δ) :
                     (kiP[busidx]+=Δ;kiM[busidx]-=Δ)
        plus=model_dc(ctx,rho,kpP,kiP);minus=model_dc(ctx,rho,kpM,kiM)
        dh=(plus.h0.-minus.h0)./(2Δ)
        da=(plus.alpha-minus.alpha)/(2Δ)
        maxabs=maximum(abs,dh);gain=="Kp" ? (maxkp=max(maxkp,maxabs)) :
                                               (maxki=max(maxki,maxabs))
        maxalpha=max(maxalpha,abs(da))
        push!(rows,(;design=j,bus=29+busidx,gain_tested=gain,rho=join(round.(rho,digits=7),";"),
            Kp=join(round.(kp,digits=5),";"),Ki=join(round.(ki,digits=5),";"),
            dH0_max_abs_Hz_per_MW_per_gain= maxabs,
            dH0_max_rel=maxabs/max(maximum(abs,base.h0),1e-12),
            dalpha_per_gain=da,alpha=base.alpha,
            centered_step=Δ,conditioning_A=cond(base.red.A)))
        println("F3_DC_POINT=",j," gain=",gain," max_abs_dH0=",maxabs)
    end
    CSV.write(joinpath(OUT,"TABLE_Q2_F3_PLL_DC_support.csv"),DataFrame(rows))
    # Dynamic derivatives on the interior corrected ExpG point. Gains are
    # moved only for diagnosis, inside their frozen ExpK box.
    cg=TOML.parsefile(joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml"))
    grow=sort!(cg["generator"],by=x->Int(x["bus"]))
    rho=Float64.([x["rho"] for x in grow]);kp=Float64.([x["Kp"] for x in grow]);ki=Float64.([x["Ki"] for x in grow])
    kidx=1;hkp=1e-3*(ctx.kpmax[kidx]-ctx.kpmin[kidx]);hki=1e-3*(ctx.kimax[kidx]-ctx.kimin[kidx])
    t=collect(0.0:0.01:8.0)
    omega_peak=Float64(ImportCSV_G07())
    function dynamic_metrics(kpp,kii)
        dm=model_dc(ctx,rho,kpp,kii)
        sig=GridFrequency.grid_step_signals(ctx,dm.m,rho,t;load_bus=16,
            disturbance_MW=1.0,gauge_vector=N.gauge_vector)
        met=GridFrequency.grid_step_metrics(sig;half_windows=(5,),dt_s=0.01)[1]
        rb=beta_at_frequency(dm.red.A,omega_peak)
        (;alpha=dm.alpha,beta_star_at_frozen_peak=rb.beta,F_peak=met.F_peak_Hz,
          R_peak=met.R_peak_Hz_s,F_inf=maximum(abs,dm.h0))
    end
    kpP=copy(kp);kpM=copy(kp);kpP[kidx]+=hkp;kpM[kidx]-=hkp
    kiP=copy(ki);kiM=copy(ki);kiP[kidx]+=hki;kiM[kidx]-=hki
    p0=dynamic_metrics(kp,ki);pp=dynamic_metrics(kpP,ki);pm=dynamic_metrics(kpM,ki)
    ip=dynamic_metrics(kp,kiP);im=dynamic_metrics(kp,kiM)
    dynrows=NamedTuple[]
    for g in ("Kp","Ki")
        h=g=="Kp" ? hkp : hki;a=g=="Kp" ? pp : ip;b=g=="Kp" ? pm : im
        for metric in (:alpha,:beta_star_at_frozen_peak,:F_peak,:R_peak,:F_inf)
            derivative=(getproperty(a,metric)-getproperty(b,metric))/(2h)
            push!(dynrows,(;point="corrected_ExpG_interior",bus=30,gain=g,metric=String(metric),
                baseline=getproperty(p0,metric),centered_derivative=derivative,step=h,
                kp=kp[1],ki=ki[1],omega_rad_s=omega_peak,
                beta_method="pointwise 2-norm radius at the frozen ExpG sampled peak; local envelope diagnostic"))
        end
    end
    CSV.write(joinpath(OUT,"TABLE_Q2_F3_dynamic_gain_derivatives.csv"),DataFrame(dynrows))
    status=maxkp<=1e-8 && maxki<=1e-8 ? "PLL_DC_SUPPORT_INDEPENDENT" : "PLL_DC_SUPPORT_PRESENT"
    result=Dict("stage"=>"F3","withheld_mixed_points"=>30,"centered_DC_gain_derivatives"=>60,
        "max_abs_dH0_dKp"=>maxkp,"max_abs_dH0_dKi"=>maxki,"max_abs_dalpha_dK"=>maxalpha,
        "PLL_STEADY_SUPPORT_RESULT"=>status,"dynamic_gain_point"=>"ExpG corrected interior",
        "dynamic_metrics"=>"alpha, sampled beta-star, SG-filtered h=5 F_peak and R_peak; F_inf",
        "beta_method"=>"pointwise fixed-frequency singular-value radius at ExpG frozen peak, not a full re-optimized beta-star or formal interval certificate")
    open(joinpath(OUT,"F3_RESULTS.toml"),"w") do io;TOML.print(io,result);end
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
    # F3 — PLL gains and sustained power support

    A centered finite difference over 30 withheld mixed designs perturbed Kp at 15 points and Ki at 15 points, one local gain per point. DC output is exact `-C A⁻¹B + D` from the fixed-support reduced model.

    - Largest sampled `|dH0/dKp|`: $(maxkp) Hz/MW per gain unit.
    - Largest sampled `|dH0/dKi|`: $(maxki) Hz/MW per gain unit.
    - Largest sampled `|dα/dK|`: $(maxalpha) s⁻¹ per gain unit.
    - Dynamic centered differences for α, sampled β*, h=5 filtered Fpeak/Rpeak, and F∞ at an interior corrected ExpG point are in `TABLE_Q2_F3_dynamic_gain_derivatives.csv`.
    - Classification: **$(status)**. β is a sampled resolvent estimate, not a formal frequency-domain certificate.
    """)
end

function ImportCSV_G07()
    tbl=CSV.read(joinpath(ROOT,"reports","experiment_G","tables","TABLE_G07_robustness_certificate.csv"),DataFrame)
    row=only(filter(r->r.case=="ExpG_exact_coordinate_candidate",eachrow(tbl)))
    row.omega_peak
end
main()
