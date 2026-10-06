using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A freeze required")
ctx=design_context(ROOT)
e=CSV.read(joinpath(ROOT,"reports","experiment_N","validation_cases",
    "ExpE_corrected","case_definition.csv"),DataFrame)
gainsets=[("nominal",fill(PDExactDesignN.K0P,10),fill(PDExactDesignN.K0I,10)),
          ("lowKp_highKi",ctx.kpmin,ctx.kimax),
          ("ExpE_retrimmed",Float64.(e.Kp),Float64.(e.Ki))]
epsgrid=[0.0,1e-5,3e-5,1e-4,3e-4,5e-4,7e-4,1e-3,2e-3,
         5e-3,1e-2,3e-2,0.1,0.3,1.0]
rows=NamedTuple[]
for bus in 30:39
    for (gainset,kp,ki) in gainsets
        for eps in epsgrid
            rho=ones(10);rho[bus-29]=1-eps
            sp=spectrum(ctx,rho,kp,ki)
            lam=sp.lambda[argmax(real.(sp.lambda))]
            push!(rows,(;bus,gainset,retained_fraction=eps,
                retained_SG_MW=ctx.power[bus-29]*eps,
                alpha=sp.alpha,rightmost_real=real(lam),
                rightmost_imag=imag(lam),n_finite=length(sp.lambda),
                architecture=eps==0 ? "all_GFL" : eps==1 ? "SG_only_at_bus" : "mixed_at_bus",
                feasible=sp.alpha<=-0.05))
        end
        sub=filter(r->r.bus==bus && r.gainset==gainset,rows)
        best=filter(r->r.feasible,sub)
        println("SCAN bus=",bus," gains=",gainset," feasible=",length(best),
            " min_cost=",isempty(best) ? NaN : minimum(r.retained_SG_MW for r in best));flush(stdout)
    end
end
CSV.write(joinpath(ROOT,"reports","experiment_N","TABLE_N09_single_support_scan.csv"),DataFrame(rows))
