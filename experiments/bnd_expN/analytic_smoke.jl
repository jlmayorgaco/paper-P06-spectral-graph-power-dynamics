using LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN

ctx=design_context(ROOT)
k=fill(PDExactDesignN.K0P,10)
i=fill(PDExactDesignN.K0I,10)
for (label,rho) in (("all_SG",zeros(10)),("bus38_half",[b==38 ? 0.5 : 0.0 for b in 30:39]))
    op=spectrum(ctx,rho,k,i)
    println("N_ANALYTIC_SMOKE ",label," alpha=",op.alpha,
        " states=",op.model.n_dynamic," gauge_residual=",op.gauge_residual)
end
