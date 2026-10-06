using CSV, DataFrames, LinearAlgebra, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const OUT=joinpath(ROOT,"reports","experiment_Q2B","CORE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"));include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","CoreDesign.jl"))
const N=ExpP.PDExactDesignN;const CD=CoreDesign;const SIG=.05;const BREQ=1.6991206999182038e-6
d=TOML.parsefile(joinpath(OUT,"Z_Q2B_CORE_d100_restored.toml"));ctx=N.design_context(ROOT);s=Int.(d["support"]);e=Float64.(d["epsilon"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"]);y=CD.encode(ctx,e,kp,ki)
ev=CD.design_eval(ctx,s,e,kp,ki;disturbance_MW=100.0,dt_s=.05,horizon_s=30.0)
Q=ev.Q
function beta0(y,ctx,s,Q)
    e,kp,ki=CD.decode(ctx,s,y);rho=ones(10);for (j,b) in enumerate(s);rho[b-29]=1-e[j];end
    m=N.descriptor(ctx,rho,kp,ki);Aq=transpose(Q)*m.Ared*Q
    minimum(svdvals(ComplexF64.(-(Aq+SIG*I))))
end
rows=NamedTuple[]
for h in (2e-4,2e-5)
    for j in eachindex(y)
        yp=copy(y);ym=copy(y);yp[j]+=h;ym[j]-=h
        bp=beta0(yp,ctx,s,Q);bm=beta0(ym,ctx,s,Q)
        push!(rows,(;step=h,coordinate=j,kind=j<=10 ? "epsilon" : (j<=20 ? "Kp" : "Ki"),
            d_beta0=(bp-bm)/(2h),dg_beta0=-(bp-bm)/(2h*BREQ),beta_plus=bp,beta_minus=bm))
    end
end
tbl=DataFrame(rows);CSV.write(joinpath(OUT,"TABLE_Q2B_beta_gradient_steps.csv"),tbl)
for h in (2e-4,2e-5)
    q=tbl[tbl.step.==h,:]
    println("BETA_GRAD step=",h," max_abs=",maximum(abs.(q.dg_beta0)),
        " eps_dg=",join(round.(q.dg_beta0[1:10];digits=5),","),
        " kp_max=",maximum(abs.(q.dg_beta0[11:20]))," ki_max=",maximum(abs.(q.dg_beta0[21:30])))
end
println("POINT beta_sampled=",ev.beta," alpha=",ev.alpha," maxg=",maximum(ev.g))
M=ComplexF64.(-(ev.quotient+SIG*I));sv=svdvals(M)
println("SIGMA_MIN_VALUES=",sv[end-4:end]," gap=",sv[end-1]-sv[end],
    " relative_gap=",(sv[end-1]-sv[end])/max(sv[end],eps()))
