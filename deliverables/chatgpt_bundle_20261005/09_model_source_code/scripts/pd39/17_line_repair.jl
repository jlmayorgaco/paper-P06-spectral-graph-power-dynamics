using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const V8 = collect(CANDIDATE_SG_BUSES)
const LINE_START = parse(Int, get(ENV, "PD39_LINE_START", "1"))
const LINE_END = parse(Int, get(ENV, "PD39_LINE_END", "46"))
const OUT = get(ENV, "PD39_LINE_OUT", joinpath(ROOT, "results", "PD39_V8_LINE_REPAIR.csv"))
const SUMMARY_OUT = get(ENV, "PD39_LINE_SUMMARY_OUT", joinpath(ROOT, "results", "PD39_V8_LINE_REPAIR_SUMMARY.csv"))
mkpath(dirname(OUT))

function evaluate_line(line, gamma)
    rows = NamedTuple[]
    delta = zeros(Float64, 46)
    delta[line] = 1 / gamma - 1
    for scenario in uncertainty_scenarios()
        ctrl = (scenario.pll_scale - 1, scenario.filter_scale - 1,
                scenario.current_control_scale - 1)
        result = try
            run_margin_case(build_confirmatory_network(V8; controller_delta = ctrl,
                branch_delta = delta, bounds = :discovery))
        catch err
            (status = "failed", error_type = "build_exception", error_message = sprint(showerror, err),
             equilibrium_status = "equilibrium_failed", alpha = NaN, margin = NaN, stable = false,
             equilibrium_residual = NaN)
        end
        push!(rows, (line = line, gamma = gamma, scenario = scenario.id,
            status = result.status, error_type = result.error_type, error_message = result.error_message,
            alpha = result.alpha, margin = result.margin, stable = result.stable,
            robust = isfinite(result.margin) && result.margin >= 0.05,
            equilibrium_status = result.equilibrium_status,
            equilibrium_residual = result.equilibrium_residual))
    end
    return rows
end

detail = NamedTuple[]
summary = NamedTuple[]
for line in LINE_START:LINE_END
    print("line max ", line, " ... ")
    r = evaluate_line(line, 1.25)
    append!(detail, r)
    good = all(x.robust for x in r)
    gamma_req = NaN
    if good
        lo, hi = 1.0, 1.25
        for _ in 1:10
            mid = (lo + hi) / 2
            trial = evaluate_line(line, mid)
            append!(detail, trial)
            if all(x.robust for x in trial)
                hi = mid
            else
                lo = mid
            end
        end
        gamma_req = hi
    end
    push!(summary, (line = line, gamma_max = 1.25,
        max_gamma_worst_margin = minimum(x.margin for x in r),
        max_gamma_worst_alpha = maximum(x.alpha for x in r),
        robust_at_1_25 = good, minimum_gamma = gamma_req,
        statuses_at_1_25 = join(unique(x.status for x in r), "|")))
    CSV.write(OUT, DataFrame(detail))
    CSV.write(SUMMARY_OUT, DataFrame(summary))
    println(" ", good, " worst_margin=", minimum(x.margin for x in r))
end
println("wrote line details=", OUT, " rows=", length(detail))
