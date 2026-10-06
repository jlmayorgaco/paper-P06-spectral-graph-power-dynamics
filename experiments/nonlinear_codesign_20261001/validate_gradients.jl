include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, ForwardDiff, Test
const R=ReducedDAE
const CONVENTION=isempty(ARGS) ? :legacy : Symbol(ARGS[1])
BLAS.set_num_threads(1);mkpath(R.OUT)
ctx=R.N.design_context(R.ROOT)
rho=fill(.8,10);kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
m=R.model(ctx,rho,kp,ki;bus=8,delta=100.,dc_convention=CONVENTION)
x=copy(m.x0);x.+=.00001sin.(eachindex(x))
d=R.derivatives(x,m)
J=ForwardDiff.jacobian(z->begin f=similar(z);R.rhs!(f,z,m,0.);f;end,x)
@test norm(d.Fx-J,Inf)/norm(J,Inf)<1e-10
P=ForwardDiff.jacobian(vcat(rho,log.(kp),log.(ki))) do p
    mm=merge(m,(;rho=p[1:10],kp=exp.(p[11:20]),ki=exp.(p[21:30])))
    f=similar(p,length(x));R.rhs!(f,x,mm,0.);f
end
pe=norm(d.Fp-P,Inf)/norm(P,Inf)
@test pe<1e-10
println("NONLINEAR_DERIVATIVES Fx=",norm(d.Fx-J,Inf)/norm(J,Inf)," Fp=",pe);flush(stdout)
horizon=2.;dt=.01
elapsed=@elapsed trajectory=R.tangent_simulate(m;horizon,dt)
peaks=R.tangent_peaks(trajectory)
println("TANGENT runtime=",elapsed," peaks=",[(a.name,a.peak,a.bus,a.time) for a in peaks]);flush(stdout)
errors=Float64[]
# Separate rho, log-Kp and log-Ki directions. The large gain sensitivities need no hand tuning.
for ind in (2,12,22)
    delta=1e-5;p=vcat(rho,log.(kp),log.(ki));pp=copy(p);pm=copy(p);pp[ind]+=delta;pm[ind]-=delta
    mp=R.model(ctx,pp[1:10],exp.(pp[11:20]),exp.(pp[21:30]);bus=8,delta=100.,dc_convention=CONVENTION)
    mm=R.model(ctx,pm[1:10],exp.(pm[11:20]),exp.(pm[21:30]);bus=8,delta=100.,dc_convention=CONVENTION)
    a=R.tangent_peaks(R.tangent_simulate(mp;horizon,dt,sensitivities=false))
    b=R.tangent_peaks(R.tangent_simulate(mm;horizon,dt,sensitivities=false))
    for j in 1:2
        fd=(a[j].peak-b[j].peak)/(2delta);exact=peaks[j].gradient[ind]
        rel=abs(fd-exact)/max(abs(fd),abs(exact),1e-7);push!(errors,rel)
        println("GRADIENT ind=",ind," metric=",peaks[j].name," tangent=",exact," FD=",fd," error=",rel);flush(stdout)
        @test rel<.003
    end
end
name=CONVENTION==:legacy ? "gradient_audit.toml" : "gradient_audit_physical.toml"
open(joinpath(R.OUT,name),"w") do io
    TOML.print(io,Dict("state_derivative_relative_error"=>norm(d.Fx-J,Inf)/norm(J,Inf),
        "parameter_derivative_relative_error"=>pe,"max_peak_gradient_relative_error"=>maximum(errors),
        "horizon_s"=>horizon,"step_s"=>dt,"tangent_runtime_s"=>elapsed,"dc_convention"=>string(CONVENTION),
        "scope"=>"Exact derivatives within a fixed support and smooth limiter regime of a second-order discretized nonlinear trajectory; no global optimality claim."))
end
