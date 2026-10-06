using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A freeze required")
ctx=design_context(ROOT)
single=CSV.read(joinpath(ROOT,"reports","experiment_N",
    "TABLE_N10_all_KKT_candidates.csv"),DataFrame)
eligible=single[(single.KKT .== true),:]
best=minimum(Float64.(eligible.retained_SG_MW))
rows=NamedTuple[]
for i in 1:9, j in (i+1):10
    for budget_factor in (0.8,0.9,0.99,1.0)
        budget=best*budget_factor
        for t in (0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
            eps_i=budget*t/ctx.power[i]
            eps_j=budget*(1-t)/ctx.power[j]
            rho=ones(10);rho[i]=1-eps_i;rho[j]=1-eps_j
            sp=spectrum(ctx,rho,ctx.kpmin,ctx.kimax)
            push!(rows,(;bus_i=i+29,bus_j=j+29,budget_factor,
                retained_SG_MW=budget,split_i=t,epsilon_i=eps_i,
                epsilon_j=eps_j,alpha=sp.alpha,
                feasible=sp.alpha<=-0.05,
                rightmost_imag=imag(sp.lambda[argmax(real.(sp.lambda))])))
        end
    end
    pair=filter(r->r.bus_i==i+29 && r.bus_j==j+29,rows)
    feas=filter(r->r.feasible && r.retained_SG_MW<best,pair)
    println("PAIR ",i+29,"_",j+29," better_feasible=",length(feas),
        " best_alpha=",minimum(r.alpha for r in pair));flush(stdout)
end
CSV.write(joinpath(ROOT,"reports","experiment_N",
    "TABLE_N09_pair_budget_audit.csv"),DataFrame(rows))
println("PAIR_AUDIT_DONE rows=",length(rows)," better_feasible=",
    count(r->r.feasible && r.retained_SG_MW<best,rows));flush(stdout)
