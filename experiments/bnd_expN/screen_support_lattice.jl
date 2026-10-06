using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A freeze required")
ctx=design_context(ROOT)
kkt=CSV.read(joinpath(ROOT,"reports","experiment_N","TABLE_N10_all_KKT_candidates.csv"),DataFrame)
best=minimum(Float64.(kkt.retained_SG_MW[kkt.KKT .== true]))
singlecost=Dict{Int,Float64}()
for bus in 30:39
    r=kkt[(kkt.bus .== bus) .& (kkt.gainset .== "lowKp_highKi"),:]
    singlecost[bus]=nrow(r)==0 ? Inf : minimum(Float64.(r.retained_SG_MW))
end
budget=best*(1-1e-4)
rows=NamedTuple[]
for mask in 0:1023
    buses=[i+29 for i in 1:10 if (mask >> (i-1)) & 1 == 1]
    for allocation in ("equal_MW","concentrated_MW")
        rho=ones(10)
        if !isempty(buses)
            if allocation=="equal_MW" || length(buses)==1
                for bus in buses
                    rho[bus-29]=1-budget/length(buses)/ctx.power[bus-29]
                end
            else
                leader=buses[argmin([singlecost[b] for b in buses])]
                for bus in buses
                    mw=bus==leader ? 0.9budget : 0.1budget/(length(buses)-1)
                    rho[bus-29]=1-mw/ctx.power[bus-29]
                end
            end
        end
        try
            sp=spectrum(ctx,rho,ctx.kpmin,ctx.kimax)
            push!(rows,(;mask,support=join(buses,";"),support_size=length(buses),
                allocation,retained_SG_MW=sum(ctx.power.*(1 .- rho)),
                alpha=sp.alpha,feasible=sp.alpha<=-0.05,
                n_finite=length(sp.lambda),status="OK",error=""))
        catch e
            push!(rows,(;mask,support=join(buses,";"),support_size=length(buses),
                allocation,retained_SG_MW=sum(ctx.power.*(1 .- rho)),
                alpha=NaN,feasible=false,n_finite=0,status="MODEL_EVALUATION_ERROR",
                error=sprint(showerror,e)))
        end
    end
    mask%64==63 && (println("LATTICE support_mask=",mask,
        " better_feasible=",count(r->r.feasible,rows));flush(stdout))
end
CSV.write(joinpath(ROOT,"reports","experiment_N",
    "TABLE_N09_support_lattice_budget_screen.csv"),DataFrame(rows))
println("LATTICE_DONE architectures=1024 samples=",length(rows),
    " better_feasible=",count(r->r.feasible,rows),
    " errors=",count(r->r.status!="OK",rows));flush(stdout)
