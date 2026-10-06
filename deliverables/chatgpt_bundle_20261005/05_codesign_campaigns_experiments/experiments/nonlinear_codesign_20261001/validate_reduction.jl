include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, ForwardDiff, Test
const R=ReducedDAE
BLAS.set_num_threads(1)
mkpath(R.OUT)
ctx=R.N.design_context(R.ROOT)
d=TOML.parsefile(joinpath(R.ROOT,"reports","codesign_validation_20261001","candidate_improved.toml"))
m0=R.model(ctx,d["rho"],d["Kp"],d["Ki"])
dx=zeros(length(m0.x0));R.rhs!(dx,m0.x0,m0,0.)
J=ForwardDiff.jacobian(x->begin y=similar(x);R.rhs!(y,x,m0,0.);y;end,m0.x0)
an=R.N.descriptor(ctx,d["rho"],d["Kp"],d["Ki"])
err=norm(J-an.Ared,Inf)/norm(an.Ared,Inf)
@test norm(dx,Inf)<1e-8
@test err<1e-10
println("EXACT_REDUCTION trim=",norm(dx,Inf)," Jacobian_relative_error=",err," n=",length(dx));flush(stdout)
rows=NamedTuple[]
for case in [(label="improved_prefix8",rho=d["rho"],kp=d["Kp"],ki=d["Ki"],bus=8,delta=100.,horizon=5.),
             (label="uniform80_bus8",rho=fill(.8,10),kp=fill(R.N.K0P,10),ki=fill(R.N.K0I,10),bus=8,delta=100.,horizon=60.)]
    m=R.model(ctx,case.rho,case.kp,case.ki;bus=case.bus,delta=case.delta)
    println("START ",case.label," vset2=",m.net.vset2);flush(stdout)
    r=R.simulate(m;horizon=case.horizon,dt=.005,tol=1e-9,wall_limit=120.)
    met=R.metrics(m,r;dt=.005,savepath=joinpath(R.OUT,case.label*".csv"))
    row=(;label=case.label,met...);push!(rows,row)
    CSV.write(joinpath(R.OUT,"reduction_cases.csv"),DataFrame(rows))
    println("RESULT ",row);flush(stdout)
end
open(joinpath(R.OUT,"reduction_identity.toml"),"w") do io
    TOML.print(io,Dict("trim_residual"=>norm(dx,Inf),"Jacobian_relative_error"=>err,
        "states"=>length(dx),"scope"=>"Nonlinear exact algebraic and passive-bus elimination; local Jacobian identity alone is not transient validation"))
end
