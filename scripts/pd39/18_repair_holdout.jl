using CSV
using DataFrames
using Statistics
using PowerDynamics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
const V8 = collect(CANDIDATE_SG_BUSES)

portfolio_from_string(s) = s == "none" ? Int[] : parse.(Int, split(s, ";"))

function voltage_range(s)
    v = Float64[]
    for bus in 1:39
        try
            push!(v, hypot(Float64(s[PowerDynamics.VIndex(bus, :busbar₊u_r)]),
                           Float64(s[PowerDynamics.VIndex(bus, :busbar₊u_i)])))
        catch
        end
    end
    isempty(v) ? (NaN, NaN) : (minimum(v), maximum(v))
end

function row_delta(r)
    [Float64(r[Symbol("delta_line_$(lpad(l, 2, '0'))")]) for l in 1:46]
end

conditions = CSV.read(joinpath(RESULTS, "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
lines = CSV.read(joinpath(RESULTS, "PD39_V8_LINE_REPAIR_SUMMARY.csv"), DataFrame)
lines = filter(r -> lowercase(string(r.robust_at_1_25)) == "true" && isfinite(Float64(r.minimum_gamma)), lines)

detail = NamedTuple[]
summary = NamedTuple[]
for lr in eachrow(lines)
    line = Int(lr.line)
    gamma = Float64(lr.minimum_gamma)
    alphas = Float64[]
    for cr in eachrow(conditions)
        ctrl = (Float64(cr.delta_pll), Float64(cr.delta_xf), Float64(cr.delta_cc))
        op = (Float64(cr.delta_load_p), Float64(cr.delta_load_q))
        # Compose the holdout line perturbation with the fixed reinforcement
        # in impedance space: Z_design = Z_base/gamma, then Z_holdout =
        # Z_design*(1+delta_holdout).
        branch_delta = row_delta(cr)
        branch_delta[line] = (1 + branch_delta[line]) / gamma - 1
        nw = build_confirmatory_network(V8; controller_delta = ctrl, load_delta = op,
            ibr_delta = Float64(cr.delta_ibr_p), branch_delta = branch_delta, bounds = :repair)
        result = try
            run_margin_case(nw)
        catch err
            (status = "failed", error_type = "exception", error_message = sprint(showerror, err),
             equilibrium_status = "equilibrium_failed", alpha = NaN, margin = NaN, stable = false,
             equilibrium_residual = NaN)
        end
        vmin, vmax = NaN, NaN
        if result.status == "ok"
            # Re-solve only the equilibrium to record voltage observables;
            # the spectrum was already evaluated by run_margin_case.
            eq = initialize_equilibrium(nw; sparse = false)
            if eq.state_finite && eq.fixed_point
                vmin, vmax = voltage_range(eq.state)
            end
        end
        isfinite(result.alpha) && push!(alphas, result.alpha)
        push!(detail, (solution = "line_$(line)", line = line, gamma = gamma,
            condition = cr.condition, status = result.status, error_type = result.error_type,
            error_message = result.error_message, alpha = result.alpha, margin = result.margin,
            stable = result.stable, robust = isfinite(result.margin) && result.margin >= 0.05,
            equilibrium_residual = result.equilibrium_residual, voltage_min_pu = vmin,
            voltage_max_pu = vmax, connected = result.status != "graph_disconnected",
            load_served = result.status == "ok"))
    end
    q = filter(r -> r.solution == "line_$(line)", detail)
    push!(summary, (solution = "line_$(line)", line = line, gamma = gamma,
        conditions = length(q), equilibrium_success = sum([r.status == "ok" for r in q]),
        robust_cases = sum([r.robust for r in q]), robust_coverage = mean([r.robust for r in q]),
        worst_alpha = isempty(alphas) ? NaN : maximum(alphas),
        median_alpha = isempty(alphas) ? NaN : median(alphas),
        worst_margin = isempty(alphas) ? NaN : -maximum(alphas),
        voltage_feasible = all([(r.voltage_min_pu >= 0.90 && r.voltage_max_pu <= 1.10) for r in q]),
        load_served = all([r.load_served for r in q])))
end

CSV.write(joinpath(RESULTS, "PD39_V8_LINE_REPAIR_HOLDOUT.csv"), DataFrame(detail))
CSV.write(joinpath(RESULTS, "PD39_V8_LINE_REPAIR_HOLDOUT_SUMMARY.csv"), DataFrame(summary))
println("line repair holdout complete")
