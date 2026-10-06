using LinearAlgebra, CSV, DataFrames, Random, Statistics
include(joinpath(@__DIR__, "..", "nonlinear_codesign_20261001", "ReducedDAE.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
const R = ReducedDAE
const OUT = @__DIR__
BLAS.set_num_threads(1)

ctx = R.N.design_context(R.ROOT)
net = PD39.ieee39_data()
theta0 = angle.(ctx.net.voltage)
vm = abs.(ctx.net.voltage)
n = length(theta0)
branches = net.branch
B = zeros(Float64, n, nrow(branches))
gamma = zeros(Float64, nrow(branches))
for (e, br) in enumerate(eachrow(branches))
    i, j = Int(br.src_bus), Int(br.dst_bus)
    y = inv(complex(Float64(br.R), Float64(br.X)))
    b = imag(-Float64(br.r_src) * y)
    B[i,e] = 1; B[j,e] = -1
    gamma[e] = b * vm[i] * vm[j]
end
eta0 = B' * theta0
weights1 = gamma .* cos.(eta0)
weights2 = -gamma .* sin.(eta0)
weights3 = -gamma .* cos.(eta0)

function branch_power(theta)
    B * (gamma .* sin.(B' * theta))
end
function taylor_terms(q)
    dq = B' * q
    p0 = branch_power(theta0)
    p1 = B * (weights1 .* dq)
    p2 = B * (weights2 .* dq.^2) / 2
    p3 = B * (weights3 .* dq.^3) / 6
    p0, p1, p2, p3
end

# Frozen validation design: 12 seeded, zero-mean directions and five amplitudes.
rng = MersenneTwister(20261002)
directions = [begin q = randn(rng,n); q .-= mean(q); q ./= norm(q); q end for _ in 1:12]
amps = [1e-4, 1e-3, 1e-2, 5e-2, 1e-1]
rows = NamedTuple[]
for (d,q0) in enumerate(directions), a in amps
    q = a .* q0
    p0,p1,p2,p3 = taylor_terms(q)
    exact = branch_power(theta0+q) - p0
    approx1 = p1
    approx2 = p1+p2
    approx3 = p1+p2+p3
    denom = max(norm(exact), eps())
    push!(rows,(;direction_id=d,amplitude_rad=a,
        first_order_relative_error=norm(exact-approx1)/denom,
        second_order_relative_error=norm(exact-approx2)/denom,
        third_order_relative_error=norm(exact-approx3)/denom,
        exact_increment_norm=norm(exact),third_order_remainder_norm=norm(exact-approx3),
        status="NUMERICALLY_VALIDATED_LOSSLESS_BRANCH_TAYLOR"))
end
CSV.write(joinpath(OUT,"TABLE_D03_TAYLOR_NETWORK_VALIDATION.csv"),DataFrame(rows))
open(joinpath(OUT,"D03_GATE_REPORT.md"),"w") do io
    println(io,"# D3 lossless network Taylor validation\n")
    println(io,"Status: NUMERICALLY_VALIDATED for the declared lossless sinusoidal branch-power map only. It does not validate lossy AC Taylor tensors, a nonlinear stability proof, or DDE frontier predictions.\n")
    println(io,"Frozen protocol: 12 seeded zero-mean angle directions, amplitudes ",join(amps,", ")," rad; 60 comparisons total. See TABLE_D03_TAYLOR_NETWORK_VALIDATION.csv.\n")
    for (name,key) in (("first_order_relative_error","first_order_relative_error"),("second_order_relative_error","second_order_relative_error"),("third_order_relative_error","third_order_relative_error"))
        vals = getproperty(DataFrame(rows), Symbol(key))
        println(io,"- ",name,": max=",maximum(vals),", median=",median(vals),".")
    end
end
println("D03 rows=",length(rows)," max third-order remainder=",maximum(r.third_order_remainder_norm for r in rows))
