using CSV, DataFrames, LinearAlgebra, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_Q2B", "CORE")
include(joinpath(ROOT, "src", "bnd_design_p", "ExpP.jl"))
include(joinpath(ROOT, "src", "bnd_expQ", "LinearSecurity.jl"))
include(joinpath(ROOT, "src", "bnd_expQ2B", "FiniteWindow.jl"))
include(joinpath(ROOT, "src", "bnd_expQ2B", "CoreDesign.jl"))
include(joinpath(ROOT, "src", "bnd_design_g", "Robustness.jl"))
const N = ExpP.PDExactDesignN

path = joinpath(OUT, "Z_Q2B_CORE_d100_restored.toml")
d = TOML.parsefile(path)
ctx = N.design_context(ROOT)
support = Int.(d["support"])
epsilon = Float64.(d["epsilon"])
rho = 1 .- epsilon
kp = Float64.(d["Kp"])
ki = Float64.(d["Ki"])
sp = N.spectrum(ctx, rho, kp, ki)
Aq = transpose(sp.quotient) * sp.model.Ared * sp.quotient
rp = resolvent_peak(Aq, 0.05, 1.0)
println("ROBUST_FULL status=", rp.status, " beta_star=", rp.beta_star,
    " peak=", rp.peak, " omega=", rp.omega_peak,
    " refinements=", rp.refinements, " grid_n=", length(rp.frequency))

metrics = FiniteWindow.design_metrics(ctx, sp.model, rho; load_bus=16,
    disturbance_MW=100.0, windows=(0.2, 0.5, 1.0, 2.0), dt_s=0.01,
    horizon_s=60.0, gauge_vector=N.gauge_vector)
CSV.write(joinpath(OUT, "TABLE_Q2B_core_rocof_windows_restored.csv"), metrics.metrics)
for r in metrics.metrics
    println("WINDOW T=", r.window_s, " Fpeak=", r.F_peak_Hz,
        " Rpeak=", r.R_peak_Hz_s, " F_inf=", r.F_inf_Hz)
end

residual = 0.0
audit = Dict{String,Any}(
    "candidate" => basename(path), "support" => support,
    "retained_SG_MW" => Float64(d["retained_SG_MW"]),
    "alpha" => sp.alpha, "beta_sampled" => Float64(d["beta_sampled"]),
    "beta_full_refined_observed" => rp.beta_star, "omega_peak" => rp.omega_peak,
    "beta_pass_observed" => rp.beta_star >= 1.6991206999182038e-6,
    "frequency_peak_T05_Hz" => only(filter(x -> x.window_s == 0.5, metrics.metrics)).F_peak_Hz,
    "rocof_peak_T05_Hz_s" => only(filter(x -> x.window_s == 0.5, metrics.metrics)).R_peak_Hz_s,
    "max_constraint_coarse" => Float64(d["max_constraint"]),
    "note" => "Robustness is deterministic full-frequency refinement, not outward-rounded formal proof; high-resolution analytic trajectory only.")
open(joinpath(OUT, "CORE_RESTORED_AUDIT.toml"), "w") do io
    TOML.print(io, audit)
end
