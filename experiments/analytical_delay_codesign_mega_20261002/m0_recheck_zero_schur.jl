using LinearAlgebra, CSV, DataFrames
include(joinpath(@__DIR__, "..", "nonlinear_codesign_20261001", "ReducedDAE.jl"))
const R = ReducedDAE
BLAS.set_num_threads(1)

ctx = R.N.design_context(R.ROOT)
rho = fill(0.875, 10)
kp = fill(R.N.K0P, 10)
ki = fill(R.N.K0I, 10)
m = R.model(ctx, rho, kp, ki; bus=8, delta=0., dc_convention=:physical_supply)
A = R.derivatives(m.x0, m).Fx
n = size(A, 1)

g, _ = R.rotation_generator(m.x0, m; jacobian=false)
g ./= norm(g)
gauge_residual = norm(A*g) / max(norm(A)*norm(g), eps())

angle_indices = Int[last(m.sgidx[i]) for i in 1:10]
E = Matrix{Float64}(I, n, n)[:, angle_indices]
Qret = Matrix(qr(E-g*(g'*E)).Q)[:, 1:10]
Qhidden = nullspace(vcat(g', Qret'))
V = hcat(Qret, Qhidden)
size(V, 2) == n-1 || error("rotational quotient rank mismatch")
Aq = V' * A * V
hidden = -Aq[11:end, 11:end]
hidden_condition = cond(hidden)
status = isfinite(hidden_condition) && hidden_condition < 1e14 ? "WELL_CONDITIONED" : "BLOCKED_SINGULAR_ZERO_FREQUENCY_SCHUR"

row = (; dynamic_dimension=n, quotient_dimension=size(Aq,1), retained_dimension=10,
    hidden_dimension=size(hidden,1), gauge_residual, hidden_condition, status)
CSV.write(joinpath(@__DIR__, "M0_ZERO_FREQUENCY_SCHUR_REPRODUCTION.csv"), DataFrame([row]))
println(row)
