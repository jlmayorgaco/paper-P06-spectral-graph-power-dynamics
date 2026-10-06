using LinearAlgebra, CSV, DataFrames, TOML, SHA, Statistics
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK
const OUT=joinpath(ROOT,"reports","experiment_K")
const TABLES=joinpath(OUT,"tables")
path=joinpath(OUT,"Z_K_NOMINAL_FINAL.toml")
strip(read(path*".sha256",String))==bytes2hex(sha256(read(path))) || error("freeze mismatch")
c=TOML.parsefile(path)
ctx=design_context(ROOT)
ev=evaluate(ctx,Float64.(c["epsilon"]),Float64.(c["Kp"]),Float64.(c["Ki"]);vectors=true)
pole=ev.lambda[first(ev.active)]
pc=BNDDesignK.BNDDesignG.ClosureSpectrum.port_closure(ctx.net,pole,
    1 .- ev.epsilon,ev.kp,ev.ki)
C=pc.C
v=svd(C).V[:,end]
k=argmax([norm(v[2i-1:2i]) for i in eachindex(ctx.buses)])
kk=collect(2k-1:2k);rr=setdiff(collect(1:20),kk)
R=C[rr,rr]\Matrix{ComplexF64}(I,18,18)
others=setdiff(collect(1:10),[k])
paths=NamedTuple[];sump=zeros(ComplexF64,2,2);sumabs=0.0
for l in others,m in others
    global sump, sumabs
    lc=collect(2l-1:2l);mc=collect(2m-1:2m)
    li=[findfirst(==(x),rr) for x in lc]
    mi=[findfirst(==(x),rr) for x in mc]
    Gklm=-C[kk,lc]*R[li,mi]*C[mc,kk]
    sump.+=Gklm;sumabs+=norm(Gklm)
    push!(paths,(pivot_bus=ctx.buses[k],from_bus=ctx.buses[l],
        via_bus=ctx.buses[m],pathway_frobenius=norm(Gklm),
        pathway_trace_real=real(tr(Gklm)),
        pathway_trace_imag=imag(tr(Gklm)),
        pole_real=real(pole),pole_imag=imag(pole)))
end
ratio=sumabs/max(norm(sump),eps())
CSV.write(joinpath(TABLES,"TABLE_K10_pairwise_pathways.csv"),DataFrame(paths))

# Y_port is an explanatory static network coordinate, not a diagonal modal
# stiffness or a substitute for the exact frequency-dependent closure.
Yreal=inv(pc.G)
Yport=zeros(ComplexF64,10,10)
for i in 1:10,j in 1:10
    Yport[i,j]=Yreal[2i-1,2j-1]+im*Yreal[2i,2j-1]
end
LG=real.(Yport);LB=imag.(Yport)
comm=LG*LB-LB*LG
commratio=norm(comm)/max(norm(LG)*norm(LB),eps())
strength=[sum(abs(LG[i,j]) for j in 1:10 if j!=i) for i in 1:10]
epsilon_strength_correlation=cor(ev.epsilon,strength)
kp_strength_correlation=std(ev.kp)>1e-12 ? cor(ev.kp,strength) : NaN
rows=NamedTuple[]
for i in 1:10
    push!(rows,(bus=ctx.buses[i],retained_fraction=ev.epsilon[i],
        Kp=ev.kp[i],Ki=ev.ki[i],LG_offdiag_strength=strength[i],
        LG_LB_commutator_ratio=commratio,
        epsilon_LG_strength_correlation=epsilon_strength_correlation,
        Kp_LG_strength_correlation=kp_strength_correlation,
        pivot_bus=ctx.buses[k],pathway_cancellation_ratio=ratio,
        interpretation="DESCRIPTIVE_GRAPH_COORDINATE_ONLY"))
end
CSV.write(joinpath(TABLES,"TABLE_K10_graph_diagnostics.csv"),DataFrame(rows))
println("BND_PIVOT_BUS: ",ctx.buses[k])
println("BND_PATHWAY_CANCELLATION_RATIO: ",ratio)
println("GRAPH_LG_LB_COMMUTATOR_RATIO: ",commratio)
