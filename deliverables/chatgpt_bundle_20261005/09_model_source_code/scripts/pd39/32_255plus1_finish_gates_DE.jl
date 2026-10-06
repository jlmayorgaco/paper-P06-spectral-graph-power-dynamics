using CSV
using DataFrames

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
const MODEL = joinpath(ROOT, "src", "pd39")
mkpath(RESULTS)

const AUDIT = CSV.read(joinpath(RESULTS, "PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv"), DataFrame)
const MODAL = CSV.read(joinpath(RESULTS, "PD39_255PLUS1_MODE_TRACKING.csv"), DataFrame)
const STATIC = CSV.read(joinpath(RESULTS, "PD39_CLASSICAL_BASELINES.csv"), DataFrame)
const DISCOVERY = CSV.read(joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)

const COMPOSITION_PAIRS = [
    (pair_id = "k6_p1", cardinality = 6, a = "30;32;33;34;36;37", b = "30;33;34;35;36;37"),
    (pair_id = "k6_p2", cardinality = 6, a = "30;32;33;34;36;38", b = "30;33;34;35;36;38"),
    (pair_id = "k6_p3", cardinality = 6, a = "30;32;33;34;37;38", b = "30;33;34;35;37;38"),
    (pair_id = "k6_p4", cardinality = 6, a = "30;32;33;36;37;38", b = "30;33;35;36;37;38"),
    (pair_id = "k6_p5", cardinality = 6, a = "30;32;34;36;37;38", b = "30;34;35;36;37;38"),
    (pair_id = "k6_p6", cardinality = 6, a = "32;33;34;36;37;38", b = "33;34;35;36;37;38"),
    (pair_id = "k7_p1", cardinality = 7, a = "32;33;34;35;36;37;38", b = "30;32;34;35;36;37;38"),
    (pair_id = "k7_p2", cardinality = 7, a = "30;32;33;34;36;37;38", b = "30;33;34;35;36;37;38"),
    (pair_id = "k7_p3", cardinality = 7, a = "30;32;33;35;36;37;38", b = "30;32;33;34;35;36;37"),
    (pair_id = "k7_p4", cardinality = 7, a = "30;32;33;34;35;37;38", b = "30;32;33;34;35;36;38"),
]

rowvalue(r, name) = getproperty(r, name)
static_row(p) = only(filter(r -> String(r.portfolio) == p, eachrow(STATIC)))
audit_row(p, s) = filter(r -> String(r.portfolio) == p && String(r.scenario) == s && String(r.audit_layer) == "independent", eachrow(AUDIT))
modal_row(p, s) = filter(r -> String(r.portfolio) == p && String(r.scenario) == s, eachrow(MODAL))
discovery_row(p, s) = filter(r -> String(r.portfolio) == p && String(r.scenario) == s, eachrow(DISCOVERY))

function composition_mechanism()
    rows = NamedTuple[]
    scenarios = unique(String.(AUDIT.scenario))
    for pair in COMPOSITION_PAIRS, sid in scenarios
        ra, rb = audit_row(pair.a, sid), audit_row(pair.b, sid)
        da, db = discovery_row(pair.a, sid), discovery_row(pair.b, sid)
        ma, mb = modal_row(pair.a, sid), modal_row(pair.b, sid)
        sa, sb = static_row(pair.a), static_row(pair.b)
        mode_ok = !isempty(ra) && !isempty(rb) && !isempty(ma) && !isempty(mb) &&
            first(ra).status == "ok" && first(rb).status == "ok" &&
            first(ma).status == "ok" && first(mb).status == "ok"
        alpha_ok = mode_ok || (!isempty(da) && !isempty(db) &&
            String(first(da).equilibrium_status) == "ok" && String(first(db).equilibrium_status) == "ok")
        push!(rows, (
            pair_id = pair.pair_id, cardinality = pair.cardinality, scenario = sid,
            portfolio_a = pair.a, portfolio_b = pair.b,
            delta_converted_mw = abs(Float64(sa.converted_mw) - Float64(sb.converted_mw)),
            delta_ibr_mva = abs(Float64(sa.ibr_mva) - Float64(sb.ibr_mva)),
            delta_remaining_sg_mw = abs(Float64(sa.remaining_sg_mw) - Float64(sb.remaining_sg_mw)),
            delta_remaining_sg_mva = abs(Float64(sa.remaining_sg_mva) - Float64(sb.remaining_sg_mva)),
            delta_remaining_inertia_mva_s = abs(Float64(sa.remaining_inertia_mva_s) - Float64(sb.remaining_inertia_mva_s)),
            alpha_a = alpha_ok ? (mode_ok ? Float64(first(ra).reference_alpha) : Float64(first(da).max_real)) : NaN,
            alpha_b = alpha_ok ? (mode_ok ? Float64(first(rb).reference_alpha) : Float64(first(db).max_real)) : NaN,
            delta_alpha_a_minus_b = alpha_ok ? ((mode_ok ? Float64(first(ra).reference_alpha) - Float64(first(rb).reference_alpha) : Float64(first(da).max_real) - Float64(first(db).max_real))) : NaN,
            mode_family_a = mode_ok ? String(first(ma).critical_mode_family) : "not_audited_6of8",
            mode_family_b = mode_ok ? String(first(mb).critical_mode_family) : "not_audited_6of8",
            mode_mac_a_to_v8 = mode_ok ? Float64(first(ma).mac_to_v8) : NaN,
            mode_mac_b_to_v8 = mode_ok ? Float64(first(mb).mac_to_v8) : NaN,
            alpha_source = mode_ok ? "AD_reference_A" : "frozen_discovery",
            status = alpha_ok ? (mode_ok ? "ok" : "ok_discovery_alpha_only") : "missing_alpha"))
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_PENETRATION_COMPOSITION.csv"), out)
    return out
end

function homotopy_audit()
    source = String[]
    for (dir, _, files) in walkdir(MODEL)
        for f in files
            endswith(f, ".jl") || continue
            push!(source, lowercase(read(joinpath(dir, f), String)))
        end
    end
    text = join(source, "\n")
    has_documented_continuous = occursin("sg_to_gfl", text) || occursin("sg-gfl interpolation", text) ||
        occursin("continuous homotopy", text)
    out = DataFrame(item = ["SG_to_GFL_physical_homotopy"],
        status = [has_documented_continuous ? "requires_manual_semantic_review" : "BLOCKED"],
        documented_continuous_path = [has_documented_continuous],
        reason = [has_documented_continuous ? "Documented continuous interpolation token found; no automatic homotopy run." :
            "Only discrete compile_bus/replace_buses semantics were found; no physically documented continuous SG-to-GFL interpolation."],
        action = ["No artificial blend constructed; no homotopy numerics run."])
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_HOMOTOPY_AUDIT.csv"), out)
    return out
end

function closure_audit()
    hits = String[]
    for (dir, _, files) in walkdir(ROOT)
        occursin(".git", dir) && continue
        for f in files
            (endswith(f, ".jl") || endswith(f, ".md") || endswith(f, ".toml")) || continue
            path = joinpath(dir, f)
            text = lowercase(read(path, String))
            if occursin("determinant identity", text) || occursin("network closure", text) ||
               occursin("q→-1", text) || occursin("q->-1", text)
                push!(hits, path)
            end
        end
    end
    evidence = isempty(hits) ?
        "No documented exact K,D,Q closure objects, dimensions, units, or determinant/Q boundary identity in the installed PD39 model/API; no proxy constructed." :
        "Repository text hits require semantic review before any exact closure claim; no proxy constructed."
    out = DataFrame(item = ["K", "D", "Q", "determinant_identity", "Q_to_minus_one_boundary"],
        status = fill("BLOCKED", 5), constructible = fill(false, 5),
        evidence = fill(evidence, 5), repository_hits = fill(join(hits, "|"), 5))
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_CLOSURE_AUDIT.csv"), out)
    return out
end

composition = composition_mechanism()
homotopy = homotopy_audit()
closure = closure_audit()
println("D composition rows=", nrow(composition), " homotopy=", homotopy.status[1])
println("E closure statuses=", unique(closure.status))
