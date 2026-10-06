using CSV, DataFrames, LinearAlgebra, SHA
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_e","PhysicalData.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveModel.jl"))

domain=PhysicalData.freeze_design_domain(ROOT)
net=CollectiveModel.frozen_network(ROOT)
buses=sort(collect(keys(net.sg))); n=length(buses)
kp=PhysicalData.nominal_pll_gains().Kp; ki=PhysicalData.nominal_pll_gains().Ki
base=CollectiveModel.spectrum(net,zeros(n),fill(kp,n),fill(ki,n))
ref=Matrix{Float64}(CSV.read(joinpath(ROOT,"reports","experiment_C","matrices","C0_Ared.csv"),DataFrame)[:,2:end])
μ=eigvals(ref)
nearest(v,w)=maximum(minimum(abs(z-y) for y in w) for z in v)
haus=max(nearest(base.lambda,μ),nearest(μ,base.lambda))
haus<1e-8 || error("frozen all-SG spectrum fails before replacement: $haus")

full=CollectiveModel.spectrum(net,ones(n),fill(kp,n),fill(ki,n))
λ=full.lambda
gauge_idx=argmin(abs.(λ))
gauge=λ[gauge_idx]
println("CLOSEST_TO_ZERO_POLE: ",gauge)
println("SMALLEST_ARED_SINGULAR_VALUE: ",minimum(svdvals(full.model.Ared)))
physical=abs(gauge)<1e-8 ? [λ[i] for i in eachindex(λ) if i!=gauge_idx] : collect(λ)
critical=physical[argmax(real.(physical))]
margin=maximum(real.(physical))
sigma=0.05
status=margin<=-sigma ? "FEASIBLE_AT_NOMINAL_GAINS" : "INFEASIBLE_AT_NOMINAL_GAINS"

dir=joinpath(ROOT,"reports","experiment_E","tables");mkpath(dir)
CSV.write(joinpath(dir,"TABLE_E04_full_replacement_nominal_spectrum.csv"),
    DataFrame(real=real.(λ),imag=imag.(λ),frequency_Hz=abs.(imag.(λ))./(2pi),
        gauge=[abs(gauge)<1e-8 && i==gauge_idx for i in eachindex(λ)],
        stable_margin=[abs(gauge)<1e-8 && i==gauge_idx ? missing : real(λ[i])<=-sigma for i in eachindex(λ)]))
CSV.write(joinpath(dir,"TABLE_E05_full_replacement_nominal_summary.csv"),DataFrame(
    domain_SHA256=[domain.sha256],generators=[n],actual_initial_SG_MW=[sum(domain.generators.SG_dispatch_initial_MW)],
    dynamic_dimension=[full.model.n_dynamic],algebraic_dimension=[full.model.n_algebraic],
    baseline_spectral_hausdorff=[haus],full_Gy_condition=[full.model.condition_Gy],
    gauge_real=[real(gauge)],gauge_imag=[imag(gauge)],
    critical_real=[real(critical)],critical_imag=[imag(critical)],
    spectral_abscissa_without_gauge=[margin],sigma_required=[sigma],
    unstable_count_without_gauge=[count(real.(physical).>0)],status=[status]))
println("DOMAIN_SHA256: ",domain.sha256)
println("INITIAL_SG_DISPATCH_MW: ",sum(domain.generators.SG_dispatch_initial_MW))
println("BASELINE_SPECTRAL_HAUSDORFF: ",haus)
println("FULL_REPLACEMENT_DYNAMIC_DIM: ",full.model.n_dynamic)
println("FULL_REPLACEMENT_GY_CONDITION: ",full.model.condition_Gy)
println("FULL_REPLACEMENT_GAUGE: ",gauge)
println("FULL_REPLACEMENT_CRITICAL_POLE: ",critical)
println("FULL_REPLACEMENT_ABSCISSA: ",margin)
println("FULL_REPLACEMENT_UNSTABLE_COUNT: ",count(real.(physical).>0))
println("FULL_REPLACEMENT_NOMINAL_STATUS: ",status)
