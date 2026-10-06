using LinearAlgebra, NetworkDynamics, PowerDynamics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
using .PDReferenceN

base=frozen_baseline()
r=zeros(10);r[31-29]=0.5
nw=build_architecture(base,r)
for rho in (0.5,0.99)
    r[31-29]=rho
    s=trim_state(nw,base,r)
    du=similar(uflat(s));nw(du,uflat(s),pflat(s),0.0)
    syms=string.(NetworkDynamics.SII.variable_symbols(nw))
    ids=sortperm(abs.(du);rev=true)[1:10]
    println("DIAG_RHO ",rho)
    for j in ids
        println(j," ",syms[j]," residual=",du[j]," state=",uflat(s)[j])
    end
end
