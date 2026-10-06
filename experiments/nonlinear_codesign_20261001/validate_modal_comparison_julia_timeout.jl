include("CoupledDQReference.jl")
using .CoupledDQReference, LinearAlgebra, DelimitedFiles, Random, TOML, SHA
using SciMLBase, OrdinaryDiffEqRosenbrock, ForwardDiff
const Q=CoupledDQReference
include("CoupledDQJacobian.jl")
BLAS.set_num_threads(1)
const OD=Q.OUT
H=readdlm(joinpath(OD,"validation_H.csv"),',',Float64)
Ti=readdlm(joinpath(OD,"validation_T.csv"),',',Float64)
center=vec(readdlm(joinpath(OD,"validation_center.csv"),',',Float64))
radii=vec(readdlm(joinpath(OD,"validation_radii.csv"),',',Float64))
group=Int.(vec(readdlm(joinpath(OD,"validation_groups.csv"),',',Float64)))
epsilon=only(readdlm(joinpath(OD,"validation_epsilon.csv"),',',Float64))
groups=[findall(==(k),group) for k in eachindex(radii)]
scale=radii[group];rng=MersenneTwister(937161)
@assert norm(H*Ti-I,Inf)<1e-8
# Integrate normalized modal coordinates, retaining all 281 physical/instrument
# states through the invertible coordinate map. No fast-state elimination.
L=H./scale;T=Ti.*scale'
raw0=Q.full_rhs(center,zeros(88)).out
println("FLOAT_REFERENCE_BASELINE metric scaled drift ",norm(L*raw0,Inf));flush(stdout)
Jschur=dq_augmented_jacobian(center,zeros(88))
Jad=ForwardDiff.jacobian(z->Q.full_rhs(z,zeros(88)).out,center)
jacerr=norm(Jschur-Jad,Inf);@assert jacerr<1e-6
println("INDEPENDENT_SCHUR_JACOBIAN max row-sum difference ",jacerr);flush(stdout)
# The small floating trim residue is not subtracted or fitted away.
function input(t,kind)
    if kind==:waves
        q=[sin((0.13+10.0^(3k/88-1))*t+0.37k)+0.3cos(0.7k+0.41t) for k in 1:88]
    else
        k=floor(Int,t/.5)
        q=[sin(1.27j+2.41k)+cos(.73j-.89k) for j in 1:88]
    end
    .8epsilon.*q./norm(q)
end
results=Dict[]
for (name,kind,boundary) in (("origin_waves",:waves,false),("origin_switches",:switches,false),("boundary_waves",:waves,true))
    y0=zeros(Q.nx)
    if boundary
        for g in groups;y0[g].=randn(rng,length(g));y0[g].*=0.85/norm(y0[g]);end
    end
    function f!(dy,y,p,t)
        z=center+T*y
        dy.=L*Q.full_rhs(z,input(t,kind)).out
    end
    function jac!(J,y,p,t)
        J.=L*dq_augmented_jacobian(center+T*y,input(t,kind))*T
    end
    start=time();watch=DiscreteCallback((u,t,integ)->time()-start>180.,integ->terminate!(integ);save_positions=(false,false))
    sol=solve(ODEProblem(ODEFunction(f!;jac=jac!),y0,(0.,20.)),Rodas5P();abstol=2e-8,reltol=2e-8,
              saveat=.01,dense=false,maxiters=200000,callback=watch,
              tstops=kind==:switches ? collect(.5:.5:20.) : Float64[])
    complete=SciMLBase.successful_retcode(sol.retcode) && sol.t[end]>=20.0 - 1e-9
    Vmax=0.;Fmax=0.;Rmax=0.;Vmin=Inf;Vmaxbus=-Inf;inputmax=0.
    for (t,y) in zip(sol.t,sol.u)
        z=center+T*y;w=input(t,kind);o=Q.full_rhs(z,w)
        Vmax=max(Vmax,maximum(norm(y[g]) for g in groups))
        Fmax=max(Fmax,maximum(abs,o.ff));Rmax=max(Rmax,maximum(abs,o.rr))
        vv=hypot.(o.v[1:2:end],o.v[2:2:end]);Vmin=min(Vmin,minimum(vv));Vmaxbus=max(Vmaxbus,maximum(vv))
        inputmax=max(inputmax,norm(w)/epsilon)
    end
    row=Dict("name"=>name,"complete"=>complete,"retcode"=>string(sol.retcode),"end_time"=>sol.t[end],
        "seconds"=>time()-start,"max_region_usage"=>Vmax,"frequency_max_Hz"=>Fmax,"rocof_max_Hz_s"=>Rmax,
        "voltage_min"=>Vmin,"voltage_max"=>Vmaxbus,"input_ball_usage"=>inputmax,
        "passed"=>complete && Vmax<1 && Fmax<.5 && Rmax<.5 && Vmin>.9 && Vmaxbus<1.1)
    push!(results,row);println(row);flush(stdout)
    open(joinpath(OD,"nonlinear_validation_julia.toml"),"w") do io
        TOML.print(io,Dict("status"=>all(r["passed"] for r in results) ? "TRAJECTORY_CHECKS_PASSED_SO_FAR" : "TRAJECTORY_CHECK_FAILED",
            "cases"=>results,"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
            "reference_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"CoupledDQReference.jl")))),
            "jacobian_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"CoupledDQJacobian.jl")))),
            "jacobian_independent_max_row_sum_error"=>jacerr,
            "certificate_sha256"=>bytes2hex(sha256(read(joinpath(OD,"block_modal_certificate.npz")))),
            "scope"=>"Three finite nonlinear trajectory checks; universal guarantee, if claimed, comes from block inequalities, not these simulations"))
    end
end
