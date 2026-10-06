using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_opt_expN","SingleSupportKKT.jl"))
using .SingleSupportKKT
const D=SingleSupportKKT.PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A freeze required")
ctx=D.design_context(ROOT)
function load_gains(name)
    t=CSV.read(joinpath(ROOT,"reports","experiment_N","validation_cases",
        name,"case_definition.csv"),DataFrame)
    Float64.(t.Kp),Float64.(t.Ki)
end
ekp,eki=load_gains("ExpE_corrected")
kkp,kki=load_gains("ExpK_corrected")
gainsets=[("nominal",fill(D.K0P,10),fill(D.K0I,10)),
          ("lowKp_highKi",ctx.kpmin,ctx.kimax),
          ("ExpE_retrimmed",ekp,eki),
          ("ExpK_retrimmed",kkp,kki)]
grid=[0.0,1e-5,3e-5,1e-4,3e-4,5e-4,7e-4,1e-3,2e-3,
      5e-3,1e-2,3e-2,0.1,0.3,1.0]
ledger=NamedTuple[];candidates=NamedTuple[]
for bus in 30:39
    for (gainset,kp,ki) in gainsets
        samples,roots=solve_single_support_branches(ctx,bus,kp,ki,grid)
        for (j,s) in enumerate(samples)
            push!(ledger,(;support_bus=bus,gainset,kind="sample",index=j,
                epsilon=s.eps,alpha=s.alpha,retained_SG_MW=ctx.power[bus-29]*s.eps,
                feasible=s.alpha<=-0.05))
        end
        for (j,r) in enumerate(roots)
            push!(ledger,(;support_bus=bus,gainset,kind="corrected_boundary",index=j,
                epsilon=r.epsilon,alpha=r.alpha,retained_SG_MW=r.retained_SG_MW,
                feasible=r.alpha<=-0.05))
            push!(candidates,(;gainset,r...))
        end
        println("BRANCH support=",bus," gains=",gainset," roots=",length(roots),
            " best_MW=",isempty(roots) ? NaN : minimum(r.retained_SG_MW for r in roots));flush(stdout)
    end
end
out=joinpath(ROOT,"reports","experiment_N")
CSV.write(joinpath(out,"TABLE_N09_architecture_support_branch_ledger.csv"),DataFrame(ledger))
CSV.write(joinpath(out,"TABLE_N10_all_KKT_candidates.csv"),DataFrame(candidates))
println("N_KKT_BRANCH_DONE samples=",length(ledger)," candidates=",length(candidates),
    " KKT=",count(r->r.KKT,candidates));flush(stdout)
