using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A freeze required")
ctx=design_context(ROOT)
rows=NamedTuple[]
for name in ("ExpE_corrected","ExpG_corrected","ExpK_corrected")
    t=CSV.read(joinpath(ROOT,"reports","experiment_N","validation_cases",name,"case_definition.csv"),DataFrame)
    rho=Float64.(t.rho);kp=Float64.(t.Kp);ki=Float64.(t.Ki)
    sp=spectrum(ctx,rho,kp,ki)
    d=simple_mode_sensitivities(ctx,rho,kp,ki)
    println(name," alpha=",sp.alpha," cost=",sum(ctx.power.*(1 .- rho)),
        " pole=",d.lambda," grad_rho=",real.(d.rho));flush(stdout)
    println(" grad_Kp=",real.(d.Kp)," grad_Ki=",real.(d.Ki));flush(stdout)
    push!(rows,(;name,alpha=sp.alpha,retained_SG_MW=sum(ctx.power.*(1 .- rho)),
        GFL_MW=sum(ctx.power.*rho),support=join((30:39)[rho .< 1],";"),
        active_real=real(d.lambda),active_imag=imag(d.lambda),
        gradient_condition=d.condition))
end
CSV.write(joinpath(ROOT,"reports","experiment_N","TABLE_N09_seed_diagnostics.csv"),DataFrame(rows))
