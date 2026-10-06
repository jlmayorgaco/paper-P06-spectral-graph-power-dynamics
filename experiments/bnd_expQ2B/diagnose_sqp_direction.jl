using LinearAlgebra, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const OUT=joinpath(ROOT,"reports","experiment_Q2B","CORE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"));include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","CoreDesign.jl"))
const N=ExpP.PDExactDesignN;const CD=CoreDesign
d=TOML.parsefile(joinpath(OUT,"Z_Q2B_CORE_d100_extended.toml"));ctx=N.design_context(ROOT);s=Int.(d["support"]);e=Float64.(d["epsilon"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"]);y=CD.encode(ctx,e,kp,ki)
ev=CD.design_eval(ctx,s,e,kp,ki;disturbance_MW=100.0,dt_s=.05,horizon_s=30.0)
jac=CD.constraint_jacobian(ctx,s,y,ev;disturbance_MW=100.0,dt_s=.05,horizon_s=30.0,fd_step=5e-4,active_tol=.5)
lo,hi=CD._bounds(ctx,s);c=zeros(length(y));for (j,b) in enumerate(s);c[j]=ctx.power[b-29]/1000;end
for trust in (.005,.001,.0002)
 A=vcat(jac.J,Matrix{Float64}(I,length(y),length(y)),-Matrix{Float64}(I,length(y),length(y)),Matrix{Float64}(I,length(y),length(y)),-Matrix{Float64}(I,length(y),length(y)))
 b=vcat(-ev.g,hi-y,y-lo,fill(trust,length(y)),fill(trust,length(y)))
 qp=CD._solve_qp(Matrix{Float64}(I,length(y),length(y)),c,A,b)
 println("DIR trust=",trust," status=",qp.status," norm=",norm(qp.d,Inf)," W=",qp.W," mu=",qp.mu," residual=",qp.residual," predg=",ev.g+jac.J*qp.d," costchange=",dot(c,qp.d))
 for a in (1.0,.5,.25,.125)
  et,kpt,kit=CD.decode(ctx,s,y+a*qp.d); evt=CD.design_eval(ctx,s,et,kpt,kit;disturbance_MW=100.0,dt_s=.05,horizon_s=30.0)
  println("   line=",a," cost=",sum(ctx.power[b-29]*et[j] for (j,b) in enumerate(s))," maxg=",maximum(evt.g)," g=",evt.g)
 end
end
