using CSV, DataFrames, LinearAlgebra, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F3")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(ROOT,"src","bnd_expQ2","GridFrequency.jl"))
include(joinpath(ROOT,"src","bnd_design_g","Robustness.jl"))
const N=ExpP.PDExactDesignN
function calc(ctx,rho,kp,ki,t,omega)
    m=N.descriptor(ctx,rho,kp,ki);red=LinearSecurity.reduced_step_model(ctx,m,rho,16;gauge_vector=N.gauge_vector)
    alpha=maximum(real.(eigvals(red.A)))
    sig=GridFrequency.grid_step_signals(ctx,m,rho,t;load_bus=16,disturbance_MW=1.0,gauge_vector=N.gauge_vector)
    met=GridFrequency.grid_step_metrics(sig;half_windows=(5,),dt_s=0.01)[1]
    R=(im*omega*I-(red.A+0.05I))\I
    beta=1/opnorm(R,2)
    (;alpha,beta,F_peak=met.F_peak_Hz,R_peak=met.R_peak_Hz_s,F_inf=maximum(abs,real.(red.C_bus*(-(red.A\red.B)).+red.D_bus)))
end
function main()
    ctx=N.design_context(ROOT);c=TOML.parsefile(joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml"));g=sort!(c["generator"],by=x->Int(x["bus"]))
    rho=Float64.([x["rho"] for x in g]);kp=Float64.([x["Kp"] for x in g]);ki=Float64.([x["Ki"] for x in g]);idx=9
    omega=CSV.read(joinpath(ROOT,"reports","experiment_G","tables","TABLE_G07_robustness_certificate.csv"),DataFrame)
    omega=only(filter(r->r.case=="ExpG_exact_coordinate_candidate",eachrow(omega))).omega_peak
    hk=Dict("Kp"=>1e-3*(ctx.kpmax[idx]-ctx.kpmin[idx]),"Ki"=>1e-3*(ctx.kimax[idx]-ctx.kimin[idx]))
    t=collect(0.0:0.01:8.0);base=calc(ctx,rho,kp,ki,t,omega);rows=NamedTuple[]
    for gain in ("Kp","Ki")
        h=hk[gain];kpP=copy(kp);kpM=copy(kp);kiP=copy(ki);kiM=copy(ki)
        if gain=="Kp";kpP[idx]+=h;kpM[idx]-=h;else;kiP[idx]+=h;kiM[idx]-=h;end
        plus=calc(ctx,rho,kpP,kiP,t,omega);minus=calc(ctx,rho,kpM,kiM,t,omega)
        for metric in (:alpha,:beta,:F_peak,:R_peak,:F_inf)
            push!(rows,(;point="ExpG_support_bus_38_GFL",bus=38,gain=gain,metric=String(metric),
                baseline=getproperty(base,metric),centered_derivative=(getproperty(plus,metric)-getproperty(minus,metric))/(2h),
                step=h,omega_rad_s=omega,beta_method="fixed ExpG active peak frequency pointwise singular-value inverse"))
        end
    end
    CSV.write(joinpath(OUT,"TABLE_Q2_F3_dynamic_gain_derivatives.csv"),DataFrame(rows))
    println("F3_EXP_G_BUS38 alpha=",base.alpha," beta_w=",base.beta," Fpeak=",base.F_peak," Rpeak=",base.R_peak)
    println("Kp derivatives=",[(r.metric,r.centered_derivative) for r in rows if r.gain=="Kp"])
    println("Ki derivatives=",[(r.metric,r.centered_derivative) for r in rows if r.gain=="Ki"])
end
main()
