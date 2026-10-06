using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
ctx=design_context(ROOT)
t=CSV.read(joinpath(ROOT,"reports","experiment_N","validation_cases","random_1","case_definition.csv"),DataFrame)
rho=Float64.(t.rho);kp=Float64.(t.Kp);ki=Float64.(t.Ki)
d=simple_mode_sensitivities(ctx,rho,kp,ki)
i=1
println("ANALYTIC ",d.rho[i]);flush(stdout)
for h in (1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,3e-6,1e-6)
    p=copy(rho);m=copy(rho);p[i]+=h;m[i]-=h
    lp=spectrum(ctx,p,kp,ki).lambda
    lm=spectrum(ctx,m,kp,ki).lambda
    vp=lp[argmin(abs.(lp .- d.lambda))]
    vm=lm[argmin(abs.(lm .- d.lambda))]
    fd=(vp-vm)/(2h)
    println("STEP ",h," fd=",fd," err=",abs(fd-d.rho[i]));flush(stdout)
end
