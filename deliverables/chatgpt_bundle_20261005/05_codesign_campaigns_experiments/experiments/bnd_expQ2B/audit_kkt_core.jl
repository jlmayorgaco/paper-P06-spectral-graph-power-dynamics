using CSV, DataFrames, LinearAlgebra, TOML
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_Q2B", "CORE")
include(joinpath(ROOT, "src", "bnd_design_p", "ExpP.jl"))
include(joinpath(ROOT, "src", "bnd_expQ", "LinearSecurity.jl"))
include(joinpath(ROOT, "src", "bnd_expQ2B", "FiniteWindow.jl"))
include(joinpath(ROOT, "src", "bnd_expQ2B", "CoreDesign.jl"))
const N = ExpP.PDExactDesignN
const CD = CoreDesign

d = TOML.parsefile(joinpath(OUT, "Z_Q2B_CORE_d100_restored.toml"))
ctx = N.design_context(ROOT)
support = Int.(d["support"])
epsilon = Float64.(d["epsilon"])
kp = Float64.(d["Kp"])
ki = Float64.(d["Ki"])
y = CD.encode(ctx, epsilon, kp, ki)
ev = CD.design_eval(ctx, support, epsilon, kp, ki;
    disturbance_MW=100.0, dt_s=0.05, horizon_s=30.0)
Jall = CD.constraint_jacobian(ctx, support, y, ev; disturbance_MW=100.0,
    fd_step=2e-4, active_tol=0.5).J
active = findall(abs.(ev.g) .<= 2e-4)
c = zeros(length(y))
for (j,b) in enumerate(support)
    c[j] = ctx.power[b-29] / 1000
end
Ja = Jall[active, :]
mu = isempty(active) ? Float64[] : pinv(transpose(Ja)) * (-c)
stationarity = isempty(active) ? norm(c, Inf) : norm(c + transpose(Ja) * mu, Inf)
evnames = ["modal", "robust_beta", "steady_frequency", "peak_frequency"]
labels = evnames[active]
data = Dict{String,Any}(
    "candidate" => basename(joinpath(OUT,"Z_Q2B_CORE_d100_extended.toml")),
    "constraint_values_scaled" => ev.g,
    "constraint_labels" => evnames,
    "near_active_indices" => active,
    "near_active_labels" => labels,
    "multipliers_unconstrained_LS" => mu,
    "stationarity_inf_norm" => stationarity,
    "dual_feasible" => all(mu .>= -1e-8),
    "LICQ_rank" => isempty(active) ? 0 : rank(Ja; atol=1e-9),
    "near_active_constraint_count" => length(active),
    "primal_violation" => max(0.0, maximum(ev.g)),
    "alpha" => ev.alpha,
    "beta_sampled" => ev.beta,
    "beta_sampled_omega_rad_s" => ev.beta_observed.omega,
    "Finf_Hz" => ev.Finf,
    "Fpeak_T05_Hz" => ev.Fpeak,
    "gradient_matrix_condition" => isempty(active) ? 1.0 : cond(Ja),
    "SOSC" => "NOT_TESTED_HESSIAN_NOT_CERTIFIED",
    "model_sha" => "e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",
)
open(joinpath(OUT,"CORE_KKT_AUDIT_RESTORED.toml"),"w") do io
    TOML.print(io,data)
end
println("KKT_VALUES=",ev.g," active=",labels," mu=",mu,
    " stat=",stationarity," rank=",data["LICQ_rank"]," cond=",data["gradient_matrix_condition"])
println("BETA_OMEGA=",ev.beta_observed.omega," beta_grad_epsilon=",Jall[2,1:10],
    " beta_grad_Kp_max=",maximum(abs.(Jall[2,11:20])),
    " beta_grad_Ki_max=",maximum(abs.(Jall[2,21:30])))
