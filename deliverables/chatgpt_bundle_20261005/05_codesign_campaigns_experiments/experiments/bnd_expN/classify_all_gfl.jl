using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("Milestone A freeze required")
ctx=design_context(ROOT)
rows=NamedTuple[]
for (name,kp,ki) in (("nominal",fill(PDExactDesignN.K0P,10),fill(PDExactDesignN.K0I,10)),
                      ("low",ctx.kpmin,ctx.kimin),("high",ctx.kpmax,ctx.kimax))
    sp=spectrum(ctx,ones(10),kp,ki)
    Aq=transpose(sp.quotient)*sp.model.Ared*sp.quotient
    F=eigen(Aq)
    j=argmax(real.(F.values));v=sp.quotient*F.vectors[:,j]
    ids=sortperm(abs.(v);rev=true)[1:5]
    states=[string(sp.model.state_inventory[i,:bus],":",
        sp.model.state_inventory[i,:kind],":",
        sp.model.state_inventory[i,:state_name]) for i in ids]
    S=svdvals(sp.model.Ared)
    Sq=svdvals(Aq)
    push!(rows,(;case="all_GFL_"*name,classification="PHYSICAL",
        subtype="GENERALIZED_ZERO_CHAIN",raw_gauge_residual=sp.gauge_residual,
        quotient_rightmost_real=real(F.values[j]),
        quotient_rightmost_imag=imag(F.values[j]),
        raw_smallest_singular=S[end],raw_second_smallest_singular=S[end-1],
        quotient_smallest_singular=Sq[end],top_states=join(states,";")))
    println(name," physical_zero=",F.values[j]," singular=",Sq[end]);flush(stdout)
end
CSV.write(joinpath(ROOT,"reports","experiment_N","TABLE_N03_all_GFL_addendum.csv"),DataFrame(rows))
