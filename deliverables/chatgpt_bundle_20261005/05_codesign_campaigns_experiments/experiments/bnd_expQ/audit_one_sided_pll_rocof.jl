using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q","Q2_REFERENCE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
using .ExpP, .LinearSecurity
const N=ExpP.PDExactDesignN
ctx=N.design_context(ROOT)
corners=[("min_min",0.25,0.25),("min_max",0.25,4.0),
         ("max_min",4.0,0.25),("max_max",4.0,4.0)]
edges=[("GFL_insertion",1e-6),("SG_removal",1-1e-6)]
rows=NamedTuple[]; rho=fill(0.0,10); kp=fill(N.K0P,10); ki=fill(N.K0I,10)
for bus in 30:39, (edge,edge_rho) in edges, (label,kpr,kir) in corners
    i=bus-29; rr=copy(rho);rr[i]=edge_rho
    pp=copy(kp);ii=copy(ki);pp[i]=kpr*N.K0P;ii[i]=kir*N.K0I
    m=N.descriptor(ctx,rr,pp,ii)
    out=frequency_outputs(m,rr)
    q=nullspace(reshape(N.gauge_vector(m)/norm(N.gauge_vector(m)),1,:))
    Bq=transpose(q)*(-m.B*(m.Gy\load_input_vector(ctx,m,16)))
    Cq=out.C*q
    ix=findfirst((Int.(out.metadata.bus).==bus).&(String.(out.metadata.kind).=="GFL"))
    ix===nothing && error("one-sided GFL output missing")
    r0=abs((Cq[ix:ix,:]*Bq)[1])
    push!(rows,(;bus,architecture_boundary=edge,gain_corner=label,Kp=pp[i],Ki=ii[i],
        rho_GFL=rr[i],epsilon_SG=1-rr[i],
        PLL_RoCoF_t0_unit_Hz_s_per_MW=r0,
        PLL_RoCoF_t0_100MW_Hz_s=100r0,
        exceeds_Rmax_100MW=100r0>0.5,
        interpretation="one-sided local boundary test at one bus; not a global support prune"))
end
CSV.write(joinpath(OUT,"TABLE_Q2_one_sided_PLL_RoCoF.csv"),DataFrame(rows))
for edge in first.(edges)
    rr=[x for x in rows if x.architecture_boundary==edge]
    println("$edge: one-sided PLL RoCoF min=$(minimum(100 .* x.PLL_RoCoF_t0_unit_Hz_s_per_MW for x in rr)) max=$(maximum(100 .* x.PLL_RoCoF_t0_unit_Hz_s_per_MW for x in rr)) over=$(count(x->x.exceeds_Rmax_100MW,rr))/$(length(rr))")
end
