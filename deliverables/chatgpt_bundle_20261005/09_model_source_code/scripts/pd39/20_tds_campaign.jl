using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
portfolio_from_string(s) = s == "none" ? Int[] : parse.(Int, split(s, ";"))

function discovery_alpha(portfolio)
    d = CSV.read(joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
    key = portfolio_string(portfolio)
    q = filter(r -> r.portfolio == key && r.scenario == "nominal", d)
    isempty(q) ? NaN : Float64(first(q).max_real)
end

function line_repair()
    p = joinpath(RESULTS, "PD39_V8_LINE_REPAIR_SUMMARY.csv")
    isfile(p) || return nothing
    d = CSV.read(p, DataFrame)
    q = filter(r -> lowercase(string(r.robust_at_1_25)) == "true" && isfinite(Float64(r.minimum_gamma)), d)
    isempty(q) ? nothing : q[argmin(Float64.(q.minimum_gamma)), :]
end

cases = NamedTuple[]
push!(cases, (case_id = "T1_best_7of8", portfolio = portfolio_from_string("30;32;33;34;35;36;38"),
    controller_delta = (0.0, 0.0, 0.0), branch_delta = zeros(Float64, 46), linear_alpha = discovery_alpha(portfolio_from_string("30;32;33;34;35;36;38"))))
push!(cases, (case_id = "T2_worst_7of8", portfolio = portfolio_from_string("30;32;33;35;36;37;38"),
    controller_delta = (0.0, 0.0, 0.0), branch_delta = zeros(Float64, 46), linear_alpha = discovery_alpha(portfolio_from_string("30;32;33;35;36;37;38"))))
push!(cases, (case_id = "T3_V8", portfolio = collect(CANDIDATE_SG_BUSES),
    controller_delta = (0.0, 0.0, 0.0), branch_delta = zeros(Float64, 46), linear_alpha = discovery_alpha(collect(CANDIDATE_SG_BUSES))))
push!(cases, (case_id = "reference_all_SG", portfolio = Int[],
    controller_delta = (0.0, 0.0, 0.0), branch_delta = zeros(Float64, 46), linear_alpha = discovery_alpha(Int[])))

lr = line_repair()
if lr !== nothing
    bd = zeros(Float64, 46)
    bd[Int(lr.line)] = 1 / Float64(lr.minimum_gamma) - 1
    push!(cases, (case_id = "T5_V8_best_line_repair", portfolio = collect(CANDIDATE_SG_BUSES),
        controller_delta = (0.0, 0.0, 0.0), branch_delta = bd,
        linear_alpha = NaN))
end

summary = NamedTuple[]
for c in cases, pulse in (0.01, 0.05)
    println(c.case_id, " pulse=", pulse)
    try
        r = simulate_tds_case(c.portfolio; pulse = pulse, controller_delta = c.controller_delta,
            branch_delta = c.branch_delta, tspan = (0.0, 20.0), saveat = 0.01)
        tag = replace(c.case_id, r"[^A-Za-z0-9_]" => "_") * "_pulse_$(Int(round(100*pulse)))pct"
        tr = DataFrame(time = r.t, voltage_pu = r.voltage,
            frequency_spread_hz = r.frequency_spread_hz,
            voltage_deviation_pu = r.voltage_deviation_pu)
        CSV.write(joinpath(RESULTS, "PD39_TDS_TRACE_$(tag).csv"), tr)
        push!(summary, (case_id = c.case_id, pulse = pulse, status = "ok",
            linear_alpha = c.linear_alpha, max_frequency_spread_hz = r.max_frequency_spread_hz,
            max_voltage_deviation_pu = r.max_voltage_deviation_pu,
            frequency_settling_s = r.frequency_settling_s,
            voltage_settling_s = r.voltage_settling_s,
            estimated_rate_s_inv = r.estimated_rate_s_inv,
            voltage_min_pu = r.voltage_min_pu, voltage_max_pu = r.voltage_max_pu,
            error = ""))
    catch err
        push!(summary, (case_id = c.case_id, pulse = pulse, status = "failed",
            linear_alpha = c.linear_alpha, max_frequency_spread_hz = NaN,
            max_voltage_deviation_pu = NaN, frequency_settling_s = NaN,
            voltage_settling_s = NaN, estimated_rate_s_inv = NaN,
            voltage_min_pu = NaN, voltage_max_pu = NaN, error = sprint(showerror, err)))
    end
end
CSV.write(joinpath(RESULTS, "PD39_TDS_SUMMARY.csv"), DataFrame(summary))
println("TDS campaign complete")
