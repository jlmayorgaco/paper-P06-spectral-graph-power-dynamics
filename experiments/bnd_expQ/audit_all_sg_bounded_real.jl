using LinearAlgebra, SHA
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q","Q2_REFERENCE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_opt_expP","RobustAudit.jl"))
using .ExpP, .RobustAudit
const N=ExpP.PDExactDesignN
req=1.6991206999182038e-6
ctx=N.design_context(ROOT)
sp=N.spectrum(ctx,zeros(10),fill(N.K0P,10),fill(N.K0I,10))
As=Matrix{Float64}(transpose(sp.quotient)*sp.model.Ared*sp.quotient)+0.05I
gamma=1/req
r=RobustAudit.hinf_upper_certificate(As,gamma)
ExpP.write_json(joinpath(OUT,"Q2_ALL_SG_BOUNDED_REAL.json"),Dict(string(k)=>v for (k,v) in pairs(r)))
println("ALLSG_CARE status=$(r.status) beta_lower=$(r.beta_lower) gamma_upper=$(r.gamma_upper) residual=$(r.care_relative_residual) Pmin=$(r.min_P_eigenvalue) LMImax=$(r.care_max_eigenvalue) closed=$(r.closed_loop_abscissa)")
