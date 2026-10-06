using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const H = CSV.read(joinpath(ROOT, "results", "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
const S = CSV.read(joinpath(ROOT, "results", "PD39_SELECTION_SET.csv"), DataFrame)
const OUT = joinpath(ROOT, "results", "PD39_HOLDOUT_CORE.csv")
mkpath(dirname(OUT))

ids = vcat(["missing_$(b)" for b in CANDIDATE_SG_BUSES],
           ["V8", "all_SG", "bus37_singleton", "representative_25pct",
            "representative_50pct", "representative_75pct"])
sel = filter(r -> r.selection_id in ids, S)
rows = NamedTuple[]
for sr in eachrow(sel)
    portfolio = sr.portfolio == "none" ? Int[] : parse.(Int, split(sr.portfolio, ";"))
    for hr in eachrow(H)
        print("holdout ", sr.selection_id, " ", hr.condition, " ... ")
        b = [Float64(hr[Symbol("delta_line_$(lpad(i, 2, '0'))")]) for i in 1:46]
        result = try
            nw = build_confirmatory_network(portfolio;
                controller_delta = (hr.delta_pll, hr.delta_xf, hr.delta_cc),
                load_delta = (hr.delta_load_p, hr.delta_load_q),
                ibr_delta = hr.delta_ibr_p, branch_delta = b)
            run_margin_case(nw)
        catch err
            (status = "failed", error_type = "build_exception", error_message = sprint(showerror, err),
             equilibrium_status = "equilibrium_failed", alpha = NaN, margin = NaN,
             stable = false, equilibrium_residual = NaN)
        end
        println(result.status, " alpha=", result.alpha)
        push!(rows, (selection_id = sr.selection_id, portfolio = sr.portfolio,
            cardinality = sr.cardinality, converted_mw = sr.converted_mw,
            condition = hr.condition, status = result.status, error_type = result.error_type,
            error_message = result.error_message, equilibrium_status = result.equilibrium_status,
            alpha = result.alpha, margin = result.margin, stable = result.stable,
            robust = isfinite(result.margin) && result.margin >= 0.05,
            equilibrium_residual = result.equilibrium_residual))
    end
    CSV.write(OUT, DataFrame(rows))
end
println("wrote ", OUT, " rows=", length(rows))
