using CSV, DataFrames, LinearAlgebra, TOML, SHA

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_Q2B", "CORE")
include(joinpath(ROOT, "src", "bnd_design_p", "ExpP.jl"))
include(joinpath(ROOT, "src", "bnd_expQ", "LinearSecurity.jl"))
include(joinpath(ROOT, "src", "bnd_expQ2B", "FiniteWindow.jl"))
include(joinpath(ROOT, "src", "bnd_expQ2B", "CoreDesign.jl"))
const N = ExpP.PDExactDesignN
const CD = CoreDesign

function write_candidate(path, cand)
    ev = cand.evaluation
    data = Dict{String,Any}(
        "stage" => "CORE_DIRECT_100MW_EXTENDED",
        "status" => cand.status,
        "support" => cand.support,
        "epsilon" => cand.epsilon,
        "rho" => cand.rho,
        "Kp" => cand.Kp,
        "Ki" => cand.Ki,
        "retained_SG_MW" => cand.retained_SG_MW,
        "converted_GFL_MW" => cand.converted_GFL_MW,
        "disturbance_MW" => 100.0,
        "alpha" => ev.alpha,
        "beta_sampled" => ev.beta,
        "beta_status" => ev.beta_observed.status,
        "F_inf_Hz" => ev.Finf,
        "F_peak_Tref_Hz" => ev.Fpeak,
        "max_constraint" => maximum(ev.g),
        "model_sha" => "e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",
    )
    open(path, "w") do io
        TOML.print(io, data)
    end
    open(path * ".sha256", "w") do io
        write(io, bytes2hex(sha256(read(path))) * "\n")
    end
end

frozen_model = TOML.parsefile(joinpath(OUT, "Z_Q2B_CORE_d100_extended.toml"))
ctx = N.design_context(ROOT)
support = Int.(frozen_model["support"])
start_epsilon = Float64.(frozen_model["epsilon"])
start_kp = Float64.(frozen_model["Kp"])
start_ki = Float64.(frozen_model["Ki"])
start_eval = CD.design_eval(ctx, support, start_epsilon, start_kp, start_ki;
    disturbance_MW=100.0, dt_s=0.05, horizon_s=30.0)
maximum(start_eval.g) <= 1e-8 || error("starting point is not feasible")
println("EXTEND_START retained=", sum(ctx.power[support .- 29] .* start_epsilon),
    " maxg=", maximum(start_eval.g), " alpha=", start_eval.alpha,
    " beta_sampled=", start_eval.beta, " Finf=", start_eval.Finf,
    " Fpeak_T05=", start_eval.Fpeak)

cand = CD.solve_fixed_support(ctx, support, start_epsilon, start_kp, start_ki;
    disturbance_MW=100.0, maxiter=50, trust_radius=0.002, dt_s=0.05,
    horizon_s=30.0, checkpoint_path=joinpath(OUT, "Z_Q2B_live_d100_restore.toml"),
    restoration=true, restoration_maxiter=4)
out_path = joinpath(OUT, "Z_Q2B_CORE_d100_restored.toml")
write_candidate(out_path, cand)
CSV.write(joinpath(OUT, "TABLE_Q2B_core_history_d100_restored.csv"), cand.history)
println("EXTEND_RESULT status=", cand.status, " retained=", cand.retained_SG_MW,
    " maxg=", maximum(cand.evaluation.g), " alpha=", cand.evaluation.alpha,
    " beta_sampled=", cand.evaluation.beta, " Finf=", cand.evaluation.Finf,
    " Fpeak_T05=", cand.evaluation.Fpeak, " accepted=", cand.accepted_steps,
    " artifact=", out_path)
