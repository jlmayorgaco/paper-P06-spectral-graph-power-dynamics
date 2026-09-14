using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const OUT = joinpath(ROOT, "results", "PD39_V8_CONTROLLER_REPAIR.csv")
mkpath(dirname(OUT))
const V8 = collect(CANDIDATE_SG_BUSES)
const GRID = collect(-0.10:0.02:0.10)

function effective(base, scale)
    return (1 + base) * scale - 1
end

function evaluate_policy(p, x, y, z)
    rows = NamedTuple[]
    for scenario in uncertainty_scenarios()
        ctrl = (effective(x, scenario.pll_scale), effective(y, scenario.filter_scale),
                effective(z, scenario.current_control_scale))
        result = try
            run_margin_case(build_confirmatory_network(V8; controller_delta = ctrl, bounds = :repair))
        catch err
            (status = "failed", error_type = "build_exception", error_message = sprint(showerror, err),
             equilibrium_status = "equilibrium_failed", alpha = NaN, margin = NaN, stable = false,
             equilibrium_residual = NaN)
        end
        push!(rows, (policy = p, delta_pll = x, delta_xf = y, delta_cc = z,
            scenario = scenario.id, status = result.status, error_type = result.error_type,
            error_message = result.error_message, alpha = result.alpha, margin = result.margin,
            stable = result.stable, robust = isfinite(result.margin) && result.margin >= 0.05,
            equilibrium_status = result.equilibrium_status,
            equilibrium_residual = result.equilibrium_residual))
    end
    return rows
end

rows = NamedTuple[]
# Full one-coordinate PLL scan is the declared first search. A 3x3 refinement
# around every PLL value is predeclared and is only used if this scan has no
# robust solution.
for x in GRID
    append!(rows, evaluate_policy("pll_only", x, 0.0, 0.0))
end
CSV.write(OUT, DataFrame(rows))

summary = combine(groupby(DataFrame(rows), [:policy, :delta_pll, :delta_xf, :delta_cc]),
    :margin => minimum => :worst_margin, :robust => all => :robust_all,
    :alpha => maximum => :worst_alpha)
if !any(summary.robust_all)
    for x in GRID, y in (-0.10, 0.0, 0.10), z in (-0.10, 0.0, 0.10)
        x == 0.0 && y == 0.0 && z == 0.0 && continue
        append!(rows, evaluate_policy("grid_3x3", x, y, z))
        CSV.write(OUT, DataFrame(rows))
    end
end

df = DataFrame(rows)
summary = combine(groupby(df, [:policy, :delta_pll, :delta_xf, :delta_cc]),
    :margin => minimum => :worst_margin, :robust => all => :robust_all,
    :alpha => maximum => :worst_alpha)
CSV.write(joinpath(ROOT, "results", "PD39_V8_CONTROLLER_REPAIR_SUMMARY.csv"), summary)
println("wrote ", OUT, " rows=", nrow(df), " policies=", nrow(summary))
