isdefined(@__MODULE__,:audit) || include("solve.jl")
using Test
@testset "Interior robust frequency and envelope derivative" begin
    for omega in (1.0,2.0,5.0)
        A=[-.2 -omega;omega -.2]
        b=L.robust_minimum(A,[0.,omega-.05])
        @test isapprox(b.beta,.15;atol=1e-11)
        @test isapprox(b.omega,omega;atol=1e-8)
        derivative=-real(dot(b.F.V[:,1],[0. -1.;1. 0.]*b.F.U[:,1]))
        @test abs(derivative)<1e-8
    end
end
d=TOML.parsefile(abspath(ARGS[1]));a=L.architecture(d["support"],d);x=L.encode(a.support,d)
v=L.evaluate(a,x;derivatives=true,fine=true);au=audit(a,x,v)
q=a.Q;g=L.N.gauge_vector(v.m);h=g/dot(g,g);xl=-q*(v.A\v.B)
bf=-v.m.B*(v.m.Gy\a.input)
dc_gauge=100dot(h,v.m.Ared*xl+bf)/(2pi)
rotation_error=norm(-a.angle*(v.m.Gy\v.m.C)*g.-ones(10),Inf)
report=Dict{String,Any}("primal"=>au.primal,"stationarity"=>au.stationarity,
    "complementarity"=>au.complementarity,"LICQ_rank"=>au.rank,"active_count"=>au.active_count,
    "active_constraints"=>au.names,"multipliers"=>au.mu,"retained_MW"=>v.J,
    "beta_frequency_stationarity"=>v.beta_frequency_stationarity,
    "alpha"=>v.alpha,"beta_observed"=>v.beta,"Fpeak"=>v.Fpeak,
    "Rpeak"=>v.Rpeak,"local_certified"=>false,
    "physical_dimension"=>size(v.A,1),"global_certified"=>false,
    "formal_certificate"=>false,"SOSC"=>"NOT_TESTED_FIRST_ORDER_GATE_FAILED",
    "eigenvector_condition"=>cond(v.V),"A_condition"=>cond(v.A),
    "dc_output_spread"=>maximum(v.dc)-minimum(v.dc),"dc_gauge_value"=>dc_gauge,
    "dc_gauge_reconstruction_error"=>maximum(abs.(v.dc.-dc_gauge)),
    "rotation_output_identity_error"=>rotation_error,
    "peak_times"=>[p.time for p in v.fp],"steady_outputs"=>v.dc,
    "stationarity_by_coordinate"=>a.c+au.A'*au.mu)
if au.primal<1e-7 && au.stationarity<1e-6 && au.complementarity<1e-7
    mu=zeros(length(v.g))
    for (label,value) in zip(au.names,au.mu)
        startswith(label,"g") && (mu[parse(Int,label[2:end])]=value)
    end
    # Keep all strictly positive-multiplier constraints as equalities. Testing
    # positivity on this larger subspace is sufficient even with weakly active
    # inequalities. Failure of this conservative test is not proof of a saddle.
    strong=findall(au.mu .>1e-8);Z=nullspace(au.A[strong,:];rtol=1e-9)
    Hessians=Matrix{Float64}[]
    for h0 in (1e-4,5e-5)
        H=zeros(length(x),length(x))
        for j in eachindex(x)
            h=min(h0,max(x[j],1-x[j])/4);xp=copy(x);xm=copy(x)
            xp[j]=min(1-1e-10,x[j]+h);xm[j]=max(1e-10,x[j]-h)
            vp=L.evaluate(a,xp;derivatives=true);vm=L.evaluate(a,xm;derivatives=true)
            H[:,j]=(vp.Jac'-vm.Jac')*mu/(xp[j]-xm[j])
        end
        push!(Hessians,Z'*((H+H')/2)*Z)
    end
    uncertainty=opnorm(Hessians[1]-Hessians[2]);dim=size(Z,2)
    mineig=dim==0 ? Inf : eigmin(Symmetric(Hessians[2]))
    report["critical_subspace_dimension"]=dim
    report["projected_Hessian_min_eigenvalue"]=mineig
    report["Hessian_step_variation_norm"]=uncertainty
    report["weakly_active_constraints"]=count(au.mu .<=1e-8)
    report["SOSC"]=(dim==0 || mineig>max(1e-8,10uncertainty)) ?
        "NUMERICALLY_POSITIVE_SUFFICIENT_SUBSPACE" : "INCONCLUSIVE_CRITICAL_CONE_OR_PRECISION"
    report["local_certified"]=false # full-band robustness and independent validation still required
end
open(joinpath(L.OUT,"SECOND_ORDER_AUDIT.toml"),"w") do io;TOML.print(io,report);end
println(report)
