include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, Test
const R=ReducedDAE
BLAS.set_num_threads(1);ctx=R.N.design_context(R.ROOT)
d=TOML.parsefile(joinpath(R.OUT,"candidate_initial_physical.toml"))
m=R.model(ctx,d["rho"],d["Kp"],d["Ki"];bus=16,delta=100.,dc_convention=:physical_supply)
x=copy(m.x0);x.+=1e-6sin.(eachindex(x));angle=.37
fx=zeros(length(x));R.rhs!(fx,x,m,0.)
y=R.rotate_state(x,m,angle);fy=zeros(length(x));R.rhs!(fy,y,m,0.)
# The derivative of a rotation shifts no angle components (the +angle is constant).
tf=R.rotate_state(fx,m,angle)
for i in 1:10;tf[last(m.sgidx[i])]-=angle;tf[m.gfidx[i][3]]-=angle;end
err=norm(fy-tf,Inf)/max(norm(fy,Inf),1.)
@test err<1e-10
rows=Any[]
for k in 1:2
    r=R.simulate(m;horizon=60.,dt=.01,tol=2e-8,wall_limit=100.,rotating_frame=true)
    met=R.metrics(m,r;dt=.01,monitor_buses=collect(1:39))
    push!(rows,met);println("ROTATING_FRAME ",met);flush(stdout)
end
open(joinpath(R.OUT,"rotating_frame_audit.toml"),"w") do io
    TOML.print(io,Dict("equivariance_relative_error"=>err,"runtime_s_second_run"=>rows[2].runtime_s,
        "Fpeak_Hz"=>rows[2].Fpeak_Hz,"Rpeak_Hz_s"=>rows[2].Rpeak_Hz_s,"complete"=>rows[2].complete,
        "scope"=>"Exact moving reference frame with phase reconstructed at all buses; no state or mode removed."))
end
