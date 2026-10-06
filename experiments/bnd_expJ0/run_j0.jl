using CSV, DataFrames, LinearAlgebra, Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "bnd_design_g", "BNDDesignG.jl"))
using .BNDDesignG

const OUT = joinpath(ROOT, "reports", "experiment_J0")
mkpath(OUT)

net = BNDDesignG.CollectiveModel.frozen_network(ROOT)
buses = sort(collect(keys(net.sg)))
n = length(buses)
nominal = BNDDesignG.PhysicalData.nominal_pll_gains()
kp = fill(nominal.Kp, n)
ki = fill(nominal.Ki, n)

function quotient(epsv)
    model = BNDDesignG.CollectiveModel.mixed_jacobian(net, 1 .- epsv, kp, ki)
    gauge = BNDDesignG.gauge_vector(net, model)
    Q = nullspace(reshape(gauge / norm(gauge), 1, :))
    Aq = transpose(Q) * model.Ared * Q
    residual = norm(model.Ared * gauge) / (norm(model.Ared) * norm(gauge))
    (; Aq, residual, values=eigvals(Aq), dynamic_dimension=model.n_dynamic)
end

base = quotient(zeros(n))
singular = sort(svdvals(base.Aq))
poles = sort(base.values; by=abs)
endpoint = DataFrame(
    dynamic_dimension=[base.dynamic_dimension],
    quotient_dimension=[size(base.Aq, 1)],
    gauge_residual=[base.residual],
    smallest_singular=[singular[1]],
    second_singular=[singular[2]],
    numerical_zero_tolerance=[1e-7],
    zeros_under_tolerance=[count(<=(1e-7), singular)],
    nearest_pole_abs=[abs(poles[1])],
    second_pole_abs=[abs(poles[2])],
    spectral_abscissa=[maximum(real.(base.values))],
)
CSV.write(joinpath(OUT, "J0_ENDPOINT.csv"), endpoint)

exponents = -8:-2
ts = 10.0 .^ collect(exponents)
function slope(x, y)
    xc = x .- mean(x)
    sum(xc .* (y .- mean(y))) / sum(abs2, xc)
end

rows = NamedTuple[]
for (i, bus) in enumerate(buses)
    roots = ComplexF64[]
    gauge_residuals = Float64[]
    dimensions = Int[]
    critical = ComplexF64[]
    for t in ts
        epsv = zeros(n)
        epsv[i] = t
        q = quotient(epsv)
        push!(roots, q.values[argmin(abs.(q.values))])
        push!(gauge_residuals, q.residual)
        push!(dimensions, q.dynamic_dimension)
        push!(critical, q.values[argmax(real.(q.values))])
    end
    logs = log10.(abs.(roots))
    q_local = slope(Float64.(collect(exponents[1:4])), logs[1:4])
    q_all = slope(Float64.(collect(exponents)), logs)
    push!(rows, (
        bus=bus, q_local_1e8_to_1e5=q_local, q_full_1e8_to_1e2=q_all,
        root_at_1e8_real=real(roots[1]), root_at_1e8_imag=imag(roots[1]),
        root_at_1e5_real=real(roots[4]), root_at_1e5_imag=imag(roots[4]),
        dimension_at_positive_epsilon=dimensions[1],
        sg_open_loop_abscissa=maximum(real.(eigvals(net.sg[bus].J.A))),
        critical_at_1e8_real=real(critical[1]),
        critical_at_1e8_imag=imag(critical[1]),
        critical_at_1e2_real=real(critical[end]),
        critical_at_1e2_imag=imag(critical[end]),
        max_gauge_residual=maximum(gauge_residuals),
    ))
end
CSV.write(joinpath(OUT, "J0_SCALING.csv"), DataFrame(rows))

println("J0 quotient zeros under 1e-7: ", endpoint.zeros_under_tolerance[1])
println("J0 local exponent range: ", extrema(getproperty.(rows, :q_local_1e8_to_1e5)))
println("J0 full exponent range: ", extrema(getproperty.(rows, :q_full_1e8_to_1e2)))
println("J0 endpoint spectral abscissa: ", endpoint.spectral_abscissa[1])
println("J0 positive-epsilon dimensions: ", unique(getproperty.(rows, :dimension_at_positive_epsilon)))
println("J0 critical-real range at 1e-8: ", extrema(getproperty.(rows, :critical_at_1e8_real)))
