using CSV, DataFrames, LinearAlgebra, TOML
include(joinpath(@__DIR__, "q3_count_dde_roots.jl"))
expdir = @__DIR__
seed = TOML.parsefile(joinpath(expdir, "seed_uniform_875.toml"))
ctx = R.N.design_context(ROOT)
L = DC.linearization(ctx, Float64.(seed["rho"]), Float64.(seed["Kp"]), Float64.(seed["Ki"]))
rows = NamedTuple[]
bal = balanced_action_model(L)
Ac, C, B = bal.Ac, bal.C, bal.B
F = schur(ComplexF64.(Ac)); T, Z = F.T, F.Z
left = C * Z; right = adjoint(Z) * B
for tau_ms in (20.0, 40.0)
    tau = fill(tau_ms / 1000, length(L.Ai))
    radius = transfer_tail_bound(Ac, B, C, tau, -0.05).radius
    step = min(20.0, pi / (4 * maximum(tau)))
    contour_points = contour(-0.05, radius, step)
    ident = Matrix{ComplexF64}(I, length(tau), length(tau))
    for (k, s) in enumerate(contour_points)
        X = (s * I - T) \ right
        EminusI = Diagonal(exp.(-s .* tau) .- 1)
        D = det(ident - EminusI * left * X)
        push!(rows, (; tau_ms, contour_index=k, s_real=real(s), s_imag=imag(s), D_real=real(D), D_imag=imag(D), D_logabs=log(abs(D)), D_phase=angle(D), tail_radius=radius))
    end
end
CSV.write(joinpath(expdir, "TABLE_Q03_REDUCED_NYQUIST_DIAGNOSTIC.csv"), DataFrame(rows))
println("saved sampled reduced Nyquist diagnostic rows=", length(rows))
