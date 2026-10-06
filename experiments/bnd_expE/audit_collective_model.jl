using CSV, DataFrames, LinearAlgebra, SHA
const ROOT = normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_e","PhysicalData.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveModel.jl"))

domain=PhysicalData.freeze_design_domain(ROOT)
println("DOMAIN_SHA256: ",domain.sha256)
net=CollectiveModel.frozen_network(ROOT)
println("NO_LOAD_KCL_INF: ",net.no_load_kcl_inf)
CSV.write(joinpath(ROOT,"reports","experiment_E","tables","TABLE_E03_initialized_ZIP_loads.csv"),net.load_audit)

n=length(net.sg); kp=PhysicalData.nominal_pll_gains().Kp; ki=PhysicalData.nominal_pll_gains().Ki
mod=CollectiveModel.mixed_jacobian(net,zeros(n),fill(kp,n),fill(ki,n))
df=CSV.read(joinpath(ROOT,"reports","experiment_C","matrices","C0_Ared.csv"),DataFrame)
aref=Matrix{Float64}(df[:,2:end]); λ=eigvals(mod.Ared); μ=eigvals(aref)
nearest(v,w)=maximum(minimum(abs(z-y) for y in w) for z in v)
haus=max(nearest(λ,μ),nearest(μ,λ))
println("BASELINE_ANALYTIC_DIM: ",mod.n_dynamic)
println("BASELINE_FROZEN_DIM: ",size(aref,1))
println("BASELINE_SPECTRAL_HAUSDORFF: ",haus)
println("BASELINE_ANALYTIC_ABSCISSA: ",maximum(real.(λ)))
println("BASELINE_FROZEN_ABSCISSA: ",maximum(real.(μ)))
println("BASELINE_ANALYTIC_LEAST_DAMPED: ",sort(λ,by=real,rev=true)[1:8])
println("BASELINE_FROZEN_LEAST_DAMPED: ",sort(μ,by=real,rev=true)[1:8])
