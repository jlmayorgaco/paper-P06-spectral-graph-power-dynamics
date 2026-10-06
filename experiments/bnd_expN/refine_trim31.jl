using LinearAlgebra, NetworkDynamics, PowerDynamics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
using .PDReferenceN

base=frozen_baseline()
r=zeros(10);r[2]=0.99
nw=build_architecture(base,r)
s=trim_state(nw,base,r)
du=similar(uflat(s));nw(du,uflat(s),pflat(s),0.0)
names=["ctrld_gen₊avr₊vm","ctrld_gen₊machine₊ψ″_q","ctrld_gen₊machine₊ψ″_d"]
idx=[findfirst(==(string(VIndex(31,Symbol(n)))),string.(NetworkDynamics.SII.variable_symbols(nw))) for n in names]
@show idx
u=copy(uflat(s));J=zeros(length(idx),length(idx));h=1e-7
for (col,j) in enumerate(idx)
    up=copy(u);um=copy(u);up[j]+=h;um[j]-=h
    fp=similar(du);fm=similar(du)
    nw(fp,up,pflat(s),0.0);nw(fm,um,pflat(s),0.0)
    J[:,col]=(fp[idx]-fm[idx])/(2h)
end
@show du[idx] J
delta= -J\du[idx]
@show delta
u[idx]+=delta
nw(du,u,pflat(s),0.0)
@show maximum(abs,du) du[idx]
for j in sortperm(abs.(du);rev=true)[1:8]
    println(j," ",NetworkDynamics.SII.variable_symbols(nw)[j]," ",du[j])
end
