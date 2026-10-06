using LinearAlgebra, TOML, Dates
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const OUT=joinpath(ROOT,"reports","experiment_Q2B","CORE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_opt_expP","RobustAudit.jl"))
const N=ExpP.PDExactDesignN;const R=RobustAudit;const BREQ=1.6991206999182038e-6
d=TOML.parsefile(joinpath(OUT,"Z_Q2B_CORE_d100_restored.toml"));ctx=N.design_context(ROOT)
s=Int.(d["support"]);eps=Float64.(d["epsilon"]);rho=1 .- eps;kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
sp=N.spectrum(ctx,rho,kp,ki);Aq=transpose(sp.quotient)*sp.model.Ared*sp.quotient;As=Matrix{Float64}(Aq)+0.05I
start=time();c=R.certify_radius(As,BREQ;max_nodes=50000,max_depth=64);elapsed=time()-start
out=Dict{String,Any}("status"=>c.status,"certified_float64_lipschitz"=>c.certified,
    "beta_lower_bound_float64"=>c.beta_lower_cert,"beta_observed_upper"=>c.beta_upper_observed,
    "beta_req"=>BREQ,"nodes"=>c.nodes,"max_depth"=>c.max_depth,
    "tail_lower"=>c.tail_lower,"elapsed_s"=>elapsed,
    "outward_rounded"=>false,"formal_certificate"=>false,
    "candidate"=>basename(joinpath(OUT,"Z_Q2B_CORE_d100_restored.toml")),
    "model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a")
open(joinpath(OUT,"CORE_BETA_INTERVAL_CERTIFICATE.toml"),"w") do io;TOML.print(io,out);end
println("BETA_INTERVAL status=",c.status," certified=",c.certified,
    " lower=",c.beta_lower_cert," observed=",c.beta_upper_observed,
    " req=",BREQ," nodes=",c.nodes," depth=",c.max_depth," elapsed_s=",elapsed)
