using CSV
using DataFrames
using Statistics

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
mkpath(RESULTS)

function qsafe(x, p)
    y = filter(isfinite, Float64.(x))
    isempty(y) ? NaN : quantile(y, p)
end

function merge_line_outputs()
    paths = filter(isfile, [joinpath(RESULTS, "PD39_V8_LINE_REPAIR_1_23.csv"),
                            joinpath(RESULTS, "PD39_V8_LINE_REPAIR_24_46.csv")])
    isempty(paths) && return
    detail = vcat([CSV.read(p, DataFrame) for p in paths]...)
    sort!(detail, [:line, :gamma, :scenario])
    CSV.write(joinpath(RESULTS, "PD39_V8_LINE_REPAIR.csv"), detail)
    spaths = filter(isfile, [joinpath(RESULTS, "PD39_V8_LINE_REPAIR_SUMMARY_1_23.csv"),
                             joinpath(RESULTS, "PD39_V8_LINE_REPAIR_SUMMARY_24_46.csv")])
    if !isempty(spaths)
        summary = vcat([CSV.read(p, DataFrame) for p in spaths]...)
        sort!(summary, :line)
        CSV.write(joinpath(RESULTS, "PD39_V8_LINE_REPAIR_SUMMARY.csv"), summary)
    end
end

function postprocess_modal()
    path = joinpath(RESULTS, "PD39_7OF8_TO_8OF8_MODAL_ANALYSIS.csv")
    isfile(path) || return
    d = CSV.read(path, DataFrame)
    rows = NamedTuple[]
    for pid in unique(d.portfolio_id)
        q = filter(r -> r.portfolio_id == pid, d)
        ok = filter(r -> r.status == "ok", q)
        isempty(ok) && continue
        worst = ok[argmax(Float64.(ok.alpha)), :]
        push!(rows, (portfolio_id = pid, portfolio = first(ok.portfolio), missing_bus = first(ok.missing_bus),
            cases = nrow(q), ok_cases = nrow(ok), worst_alpha = maximum(Float64.(ok.alpha)),
            worst_margin = minimum(Float64.(ok.margin)), critical_frequency_hz = worst.critical_frequency_hz,
            damping_ratio = worst.damping_ratio, critical_mode_family = worst.critical_mode_family,
            worst_scenario = worst.scenario, worst_eigenvector_condition = worst.eigenvector_condition,
            worst_jacobian_condition = worst.jacobian_condition,
            worst_g_z_condition = worst.g_z_condition,
            worst_smallest_singular_value = worst.smallest_singular_value,
            worst_mac_to_v8 = (pid == "V8" ? 1.0 : minimum(skipmissing(Float64.(ok.mac_to_v8)))),
            mean_mac_to_v8 = (pid == "V8" ? 1.0 : mean(skipmissing(Float64.(ok.mac_to_v8))))))
    end
    CSV.write(joinpath(RESULTS, "PD39_MODAL_SUMMARY.csv"), DataFrame(rows))
end

function postprocess_holdout()
    path = joinpath(RESULTS, "PD39_HOLDOUT_CORE.csv")
    isfile(path) || return
    d = CSV.read(path, DataFrame)
    rows = NamedTuple[]
    for sid in unique(d.selection_id)
        q = filter(r -> r.selection_id == sid, d)
        ok = filter(r -> r.status == "ok" && isfinite(r.alpha), q)
        alpha = Float64.(ok.alpha)
        margin = Float64.(ok.margin)
        push!(rows, (selection_id = sid, portfolio = first(q.portfolio), cardinality = first(q.cardinality),
            converted_mw = first(q.converted_mw), conditions = nrow(q), ok_cases = nrow(ok),
            failed_cases = nrow(q) - nrow(ok), robust_cases = sum(margin .>= 0.05),
            robust_coverage = isempty(margin) ? NaN : mean(margin .>= 0.05),
            worst_alpha = isempty(alpha) ? NaN : maximum(alpha), median_alpha = isempty(alpha) ? NaN : median(alpha),
            alpha_q1 = qsafe(alpha, .25), alpha_q3 = qsafe(alpha, .75),
            worst_margin = isempty(margin) ? NaN : minimum(margin),
            voltage_feasibility = "not_computed_in_fast_holdout"))
    end
    CSV.write(joinpath(RESULTS, "PD39_HOLDOUT_CORE_SUMMARY.csv"), DataFrame(rows))
end

function postprocess_repairs()
    rows = NamedTuple[]
    cp = joinpath(RESULTS, "PD39_V8_CONTROLLER_REPAIR_SUMMARY.csv")
    if isfile(cp)
        d = CSV.read(cp, DataFrame)
        for r in eachrow(d)
            push!(rows, (family = "controller_only", solution = r.policy,
                robust_all_discovery = r.robust_all, worst_margin = r.worst_margin,
                worst_alpha = r.worst_alpha, effort = maximum(abs.(Float64[r.delta_pll, r.delta_xf, r.delta_cc])),
                status = r.robust_all ? "pass" : "fail"))
        end
    end
    lp = joinpath(RESULTS, "PD39_V8_LINE_REPAIR_SUMMARY.csv")
    if isfile(lp)
        d = CSV.read(lp, DataFrame)
        for r in eachrow(d)
            push!(rows, (family = "line_only", solution = "line_$(r.line)",
                robust_all_discovery = r.robust_at_1_25, worst_margin = r.max_gamma_worst_margin,
                worst_alpha = r.max_gamma_worst_alpha,
                effort = Float64(r.gamma_max) - 1, status = r.robust_at_1_25 ? "pass" : "fail"))
        end
    end
    CSV.write(joinpath(RESULTS, "PD39_V8_REPAIR_SUMMARY.csv"), DataFrame(rows))
end

function exact_discovery_failures()
    path = joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv")
    isfile(path) || return
    d = CSV.read(path, DataFrame)
    q = filter(r -> r.dynamic_margin < 0.05, d)
    CSV.write(joinpath(RESULTS, "PD39_DISCOVERY_EXACT_FAILURES.csv"), q)
end

merge_line_outputs()
postprocess_modal()
postprocess_holdout()
postprocess_repairs()
exact_discovery_failures()
println("confirmatory postprocessing complete")
