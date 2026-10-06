# Direct post-freeze architecture feasibility screen, without PowerDynamics.
using CSV, DataFrames
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A must be frozen before architecture screening")
ctx=design_context(ROOT)
rows=NamedTuple[]
for (name,rho,kp,ki) in (
    ("all_GFL_nominal",ones(10),fill(PDExactDesignN.K0P,10),fill(PDExactDesignN.K0I,10)),
    ("all_GFL_low",ones(10),ctx.kpmin,ctx.kimin),
    ("all_GFL_high",ones(10),ctx.kpmax,ctx.kimax),
    ("all_SG",zeros(10),fill(PDExactDesignN.K0P,10),fill(PDExactDesignN.K0I,10)))
    sp=spectrum(ctx,rho,kp,ki)
    push!(rows,(;case=name,alpha=sp.alpha,gauge_residual=sp.gauge_residual,
        retained_SG_MW=sum(ctx.power.*(1 .- rho)),
        GFL_MW=sum(ctx.power.*rho),n_finite=length(sp.lambda),
        feasible=sp.alpha<=-0.05))
    println(name," alpha=",sp.alpha," finite=",length(sp.lambda));flush(stdout)
end
CSV.write(joinpath(ROOT,"reports","experiment_N","TABLE_N09_architecture_screen.csv"),DataFrame(rows))
