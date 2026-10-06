using CSV, DataFrames, LinearAlgebra, SHA
include(joinpath(@__DIR__, "..", "..", "src", "bnd_design_p", "ExpP.jl"))
using .ExpP
include(joinpath(@__DIR__, "..", "..", "src", "bnd_design_p", "AlgebraicFrontiers.jl"))
using .AlgebraicFrontiers

const OUT = joinpath(ExpP.ROOT, "reports", "experiment_P", "P2", "TABLE_P2_guard_policy.csv")
const SIGMA_REQ = 0.05

function main()
    ctx = ExpP.PDExactDesignN.design_context(ExpP.ROOT)
    c, candidate_sha = ExpP.read_candidate(ExpP.ROOT)
    rho = Float64.(c["rho"]); kp = Float64.(c["Kp"]); ki = Float64.(c["Ki"])
    p0 = CSV.read(joinpath(ExpP.ROOT, "reports", "experiment_P", "P0", "TABLE_P0_regression.csv"), DataFrame)
    error_est = maximum(Float64.(p0.alpha_error))
    delta_num = max(1e-6, 10error_est)
    roots = AlgebraicFrontiers.retention_roots(ctx, kp, ki, 38;
        sigma=SIGMA_REQ + delta_num, delta_num=0.0)
    isempty(roots.roots) && error("no real retention root at the proposed guarded margin")
    eps0 = 1 - rho[38-29]
    row = roots.roots[argmin(abs.(roots.roots.epsilon .- eps0)), :]
    rg = ones(10); rg[38-29] = 1 - Float64(row.epsilon)
    sp = ExpP.design_spectrum(ctx, rg, kp, ki)
    root_cost = Float64(row.retained_SG_MW)
    base_cost = dot(ctx.power, 1 .- rho)
    out = DataFrame([(; candidate_sha256=candidate_sha,
        source_candidate_unchanged=bytes2hex(sha256(read(joinpath(ExpP.ROOT, "reports", "experiment_P", "P5", "Z_P_NOMINAL_FINAL.toml")))) == strip(read(joinpath(ExpP.ROOT, "reports", "experiment_P", "P5", "Z_P_NOMINAL_FINAL.toml.sha256"), String)),
        sigma_req_per_s=SIGMA_REQ, observed_analytic_PD_alpha_error=error_est,
        delta_num_per_s=delta_num, guarded_target_per_s=SIGMA_REQ+delta_num,
        bus=38, epsilon_root=Float64(row.epsilon), rho_root=1-Float64(row.epsilon),
        retained_SG_MW= root_cost, retained_MW_increase_vs_frozen= root_cost-base_cost,
        complete_spectrum_alpha=sp.alpha, guarded_feasible=sp.alpha <= -(SIGMA_REQ+delta_num)+1e-8,
        all_physical_poles_checked=length(sp.lambda)==101,
        boundary_pole_count=count(abs.(sp.lambda .+ (SIGMA_REQ + delta_num)) .< 1e-5),
        polynomial_residual=Float64(row.polynomial_residual), TF_condition=Float64(row.TF_condition),
        candidate_was_modified=false)])
    CSV.write(OUT, out)
    println(first(eachrow(out)))
end

main()
