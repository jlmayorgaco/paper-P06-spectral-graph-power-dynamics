using CSV
using DataFrames
using LinearAlgebra
using NetworkDynamics
using PowerDynamics
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
mkpath(RESULTS)

const DISCOVERY = CSV.read(
    joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"),
    DataFrame,
)
const STATIC = CSV.read(joinpath(RESULTS, "PD39_CLASSICAL_BASELINES.csv"), DataFrame)
const SCENARIOS = uncertainty_scenarios()
const V8 = portfolio_string(CANDIDATE_SG_BUSES)
const TOTAL_CANDIDATE_MW = 4620.0
const TOTAL_CANDIDATE_MVA = 6400.0
const HIGH_PLL_SCENARIO = "pll1.2_xf0.8_cc1.2"
const TOLERANCES = [1e-8, 1e-10, 1e-12]

# This list is identical to the frozen list in
# docs/PD39_COMPOSITION_NUMERICAL_PREREG.md.  It is intentionally explicit so
# that future reruns do not reselect pairs from dynamic outcomes.
const FROZEN_PAIRS = [
    (pair_id = "k6_p1", cardinality = 6,
     a = "30;32;33;34;36;37", b = "30;33;34;35;36;37"),
    (pair_id = "k6_p2", cardinality = 6,
     a = "30;32;33;34;36;38", b = "30;33;34;35;36;38"),
    (pair_id = "k6_p3", cardinality = 6,
     a = "30;32;33;34;37;38", b = "30;33;34;35;37;38"),
    (pair_id = "k6_p4", cardinality = 6,
     a = "30;32;33;36;37;38", b = "30;33;35;36;37;38"),
    (pair_id = "k6_p5", cardinality = 6,
     a = "30;32;34;36;37;38", b = "30;34;35;36;37;38"),
    (pair_id = "k6_p6", cardinality = 6,
     a = "32;33;34;36;37;38", b = "33;34;35;36;37;38"),
    (pair_id = "k7_p1", cardinality = 7,
     a = "32;33;34;35;36;37;38", b = "30;32;34;35;36;37;38"),
    (pair_id = "k7_p2", cardinality = 7,
     a = "30;32;33;34;36;37;38", b = "30;33;34;35;36;37;38"),
    (pair_id = "k7_p3", cardinality = 7,
     a = "30;32;33;35;36;37;38", b = "30;32;33;34;35;36;37"),
    (pair_id = "k7_p4", cardinality = 7,
     a = "30;32;33;34;35;37;38", b = "30;32;33;34;35;36;38"),
]

parse_portfolio(s::AbstractString) = isempty(s) || s == "none" ? Int[] : parse.(Int, split(s, ";"))

function static_row(portfolio::AbstractString)
    rows = filter(r -> String(r.portfolio) == portfolio, eachrow(STATIC))
    isempty(rows) && error("missing static metadata for $portfolio")
    return only(rows)
end

function discovery_row(portfolio::AbstractString, scenario::AbstractString)
    rows = filter(r -> String(r.portfolio) == portfolio && String(r.scenario) == scenario,
                  eachrow(DISCOVERY))
    isempty(rows) && return nothing
    return only(rows)
end

float_or_nan(x) = x isa Missing ? NaN : Float64(x)
text_or_empty(x) = x isa Missing ? "" : String(x)

function scenario_controller(s)
    return (s.pll_scale - 1, s.filter_scale - 1, s.current_control_scale - 1)
end

function static_match(a::AbstractString, b::AbstractString)
    ra, rb = static_row(a), static_row(b)
    return (
        delta_converted_mw = abs(Float64(ra.converted_mw) - Float64(rb.converted_mw)),
        delta_ibr_mva = abs(Float64(ra.ibr_mva) - Float64(rb.ibr_mva)),
        delta_remaining_sg_mw = abs(Float64(ra.remaining_sg_mw) - Float64(rb.remaining_sg_mw)),
        delta_remaining_sg_mva = abs(Float64(ra.remaining_sg_mva) - Float64(rb.remaining_sg_mva)),
        delta_remaining_inertia_mva_s = abs(Float64(ra.remaining_inertia_mva_s) -
                                             Float64(rb.remaining_inertia_mva_s)),
        normalized_euclidean = sqrt(
            (Float64(ra.converted_mw) - Float64(rb.converted_mw))^2 / TOTAL_CANDIDATE_MW^2 +
            (Float64(ra.ibr_mva) - Float64(rb.ibr_mva))^2 / TOTAL_CANDIDATE_MVA^2 +
            (Float64(ra.remaining_inertia_mva_s) - Float64(rb.remaining_inertia_mva_s))^2 /
                61000.0^2,
        ),
    )
end

function selection_rows()
    out = NamedTuple[]
    push!(out, (pair_id = "V8", cardinality = 8, side = "V8", portfolio = V8,
                converted_mw = Float64(static_row(V8).converted_mw),
                ibr_mva = Float64(static_row(V8).ibr_mva),
                remaining_sg_mw = Float64(static_row(V8).remaining_sg_mw),
                remaining_sg_mva = Float64(static_row(V8).remaining_sg_mva),
                remaining_inertia_mva_s = Float64(static_row(V8).remaining_inertia_mva_s),
                delta_converted_mw = NaN, delta_ibr_mva = NaN,
                delta_remaining_sg_mw = NaN, delta_remaining_sg_mva = NaN,
                delta_remaining_inertia_mva_s = NaN, normalized_euclidean = NaN,
                selection_rule = "full V8"))

    for bus in CANDIDATE_SG_BUSES
        p = portfolio_string(filter(!=(bus), CANDIDATE_SG_BUSES))
        r = static_row(p)
        push!(out, (pair_id = "all_7of8", cardinality = 7, side = "missing_$(bus)", portfolio = p,
                    converted_mw = Float64(r.converted_mw), ibr_mva = Float64(r.ibr_mva),
                    remaining_sg_mw = Float64(r.remaining_sg_mw),
                    remaining_sg_mva = Float64(r.remaining_sg_mva),
                    remaining_inertia_mva_s = Float64(r.remaining_inertia_mva_s),
                    delta_converted_mw = NaN, delta_ibr_mva = NaN,
                    delta_remaining_sg_mw = NaN, delta_remaining_sg_mva = NaN,
                    delta_remaining_inertia_mva_s = NaN, normalized_euclidean = NaN,
                    selection_rule = "all eight 7-of-8 predecessors"))
    end

    for pair in FROZEN_PAIRS
        m = static_match(pair.a, pair.b)
        for (side, p) in (("A", pair.a), ("B", pair.b))
            r = static_row(p)
            push!(out, (pair_id = pair.pair_id, cardinality = pair.cardinality, side = side,
                        portfolio = p, converted_mw = Float64(r.converted_mw),
                        ibr_mva = Float64(r.ibr_mva), remaining_sg_mw = Float64(r.remaining_sg_mw),
                        remaining_sg_mva = Float64(r.remaining_sg_mva),
                        remaining_inertia_mva_s = Float64(r.remaining_inertia_mva_s),
                        delta_converted_mw = m.delta_converted_mw, delta_ibr_mva = m.delta_ibr_mva,
                        delta_remaining_sg_mw = m.delta_remaining_sg_mw,
                        delta_remaining_sg_mva = m.delta_remaining_sg_mva,
                        delta_remaining_inertia_mva_s = m.delta_remaining_inertia_mva_s,
                        normalized_euclidean = m.normalized_euclidean,
                        selection_rule = "frozen matched pair"))
        end
    end
    return DataFrame(out)
end

function reuse_selected_composition(selection)
    rows = NamedTuple[]
    portfolios = unique(String.(selection.portfolio))
    for p in portfolios, scenario in SCENARIOS
        dr = discovery_row(p, scenario.id)
        dr === nothing && error("missing frozen discovery row for $p / $(scenario.id)")
        sr = static_row(p)
        push!(rows, (
            portfolio = p, cardinality = length(parse_portfolio(p)), scenario = scenario.id,
            source = "frozen_discovery_reuse", status = text_or_empty(dr.equilibrium_status),
            error_type = text_or_empty(dr.error_type), error_message = text_or_empty(dr.error_message),
            equilibrium_status = text_or_empty(dr.equilibrium_status),
            alpha_followup = float_or_nan(dr.max_real), margin_followup = float_or_nan(dr.dynamic_margin),
            stable_followup = dr.stable isa Missing ? false : Bool(dr.stable),
            equilibrium_residual = NaN,
            alpha_discovery = float_or_nan(dr.max_real), margin_discovery = float_or_nan(dr.dynamic_margin),
            delta_alpha_followup_minus_discovery = 0.0,
            converted_mw = Float64(sr.converted_mw), ibr_mva = Float64(sr.ibr_mva),
            remaining_sg_mw = Float64(sr.remaining_sg_mw), remaining_sg_mva = Float64(sr.remaining_sg_mva),
            remaining_inertia_mva_s = Float64(sr.remaining_inertia_mva_s),
        ))
    end
    return DataFrame(rows)
end

function pair_rows(pair, composition)
    m = static_match(pair.a, pair.b)
    rows = NamedTuple[]
    for scenario in SCENARIOS
        ra = only(filter(r -> String(r.portfolio) == pair.a && String(r.scenario) == scenario.id,
                         eachrow(composition)))
        rb = only(filter(r -> String(r.portfolio) == pair.b && String(r.scenario) == scenario.id,
                         eachrow(composition)))
        ok = ra.status == "ok" && rb.status == "ok"
        push!(rows, (
            pair_id = pair.pair_id, cardinality = pair.cardinality, portfolio_a = pair.a,
            portfolio_b = pair.b, scenario = scenario.id,
            delta_converted_mw = m.delta_converted_mw, delta_ibr_mva = m.delta_ibr_mva,
            delta_remaining_sg_mw = m.delta_remaining_sg_mw,
            delta_remaining_sg_mva = m.delta_remaining_sg_mva,
            delta_remaining_inertia_mva_s = m.delta_remaining_inertia_mva_s,
            normalized_euclidean = m.normalized_euclidean,
            alpha_a = ok ? ra.alpha_followup : NaN, alpha_b = ok ? rb.alpha_followup : NaN,
            delta_alpha_a_minus_b = ok ? ra.alpha_followup - rb.alpha_followup : NaN,
            abs_delta_alpha = ok ? abs(ra.alpha_followup - rb.alpha_followup) : NaN,
            margin_a = ok ? ra.margin_followup : NaN, margin_b = ok ? rb.margin_followup : NaN,
            same_true_stability = ok ? (ra.alpha_followup < 0) == (rb.alpha_followup < 0) : false,
            same_robustness = ok ? (ra.alpha_followup <= -0.05) == (rb.alpha_followup <= -0.05) : false,
            status = ok ? "paired" : "pair_incomplete",
        ))
    end
    return DataFrame(rows)
end

function pair_summaries(pairwise)
    rows = NamedTuple[]
    for pair_id in unique(pairwise.pair_id)
        sub = filter(r -> r.pair_id == pair_id && isfinite(r.abs_delta_alpha), eachrow(pairwise))
        isempty(sub) && continue
        vals = Float64[r.abs_delta_alpha for r in sub]
        nominal = only(filter(r -> r.scenario == "nominal", eachrow(pairwise)))
        high = only(filter(r -> r.scenario == HIGH_PLL_SCENARIO, eachrow(pairwise)))
        push!(rows, (
            pair_id = pair_id, cardinality = only(sub).cardinality,
            max_abs_delta_alpha = maximum(vals), median_abs_delta_alpha = median(vals),
            nominal_delta_alpha_a_minus_b = nominal.delta_alpha_a_minus_b,
            high_pll_delta_alpha_a_minus_b = high.delta_alpha_a_minus_b,
            n_scenarios = length(vals), n_true_classification_changes = count(!, Bool[r.same_true_stability for r in sub]),
            n_robustness_classification_changes = count(!, Bool[r.same_robustness for r in sub]),
            delta_converted_mw = only(sub).delta_converted_mw,
            delta_ibr_mva = only(sub).delta_ibr_mva,
            delta_remaining_inertia_mva_s = only(sub).delta_remaining_inertia_mva_s,
            normalized_euclidean = only(sub).normalized_euclidean,
        ))
    end
    return DataFrame(rows)
end

function eq_residual(state)
    nw = extract_nw(state)
    du = zeros(Float64, length(uflat(state)))
    nw(du, uflat(state), pflat(state), state.t)
    return maximum(abs, du)
end

function reduced_spectrum(red)
    ee = eigen(red.A)
    λ = ComplexF64.(ee.values)
    keep = findall(abs.(λ) .> 1e-8)
    isempty(keep) && error("no non-gauge eigenvalues")
    k = keep[argmax(real.(λ[keep]))]
    return (alpha = real(λ[k]), critical = λ[k], count = length(λ), gauge_count = length(λ) - length(keep))
end

function initialize_with_tolerance(nw, tol)
    pf = solve_powerflow(nw; verbose = false, sparse = false,
                         tol = tol, abstol = tol, reltol = tol)
    state = initialize_from_pf!(nw; pfs = pf, verbose = false, sparsepf = false, check = :none)
    return (state = state, pf_finite = all(isfinite, uflat(pf)),
            state_finite = all(isfinite, uflat(state)),
            fixed_point = isfixpoint(state; tol = tol), residual = eq_residual(state))
end

function tolerance_audit_row(portfolio, scenario, tol)
    c = scenario_controller(scenario)
    base = (portfolio = portfolio, scenario = scenario.id, audit_type = "solver_tolerance",
            eig_method = "reference_eigen", arithmetic = "Float64", tolerance = tol,
            status = "failed", error_type = "", error_message = "", alpha = NaN,
            critical_real = NaN, critical_imag = NaN, equilibrium_residual = NaN,
            fixed_point = false, jacobian_condition = NaN, smallest_singular_value = NaN,
            generalized_alpha = NaN, descriptor_alpha_difference = NaN,
            finite_generalized_count = 0, reference_alpha = NaN, alpha_difference = NaN)
    try
        nw = build_confirmatory_network(parse_portfolio(portfolio); controller_delta = c, bounds = :discovery)
        eq = initialize_with_tolerance(nw, tol)
        sys = linearize_network(eq.state)
        red = reduce_dae(sys)
        sp = reduced_spectrum(red)
        cond = PD39._conditioning(sys, red)
        return merge(base, (status = "ok", alpha = sp.alpha, critical_real = real(sp.critical),
            critical_imag = imag(sp.critical), equilibrium_residual = eq.residual,
            fixed_point = eq.fixed_point, jacobian_condition = cond.jacobian_condition,
            smallest_singular_value = cond.smallest_singular_value, reference_alpha = sp.alpha,
            alpha_difference = 0.0))
    catch err
        return merge(base, (error_type = "exception", error_message = sprint(showerror, err)))
    end
end

function spectrum_rows(portfolio, scenario; precision = false)
    c = scenario_controller(scenario)
    rows = NamedTuple[]
    base = (portfolio = portfolio, scenario = scenario.id, audit_type = "eigensolver_crosscheck",
            tolerance = NaN, status = "failed", error_type = "", error_message = "",
            eig_method = "", arithmetic = "", alpha = NaN, critical_real = NaN,
            critical_imag = NaN, equilibrium_residual = NaN, fixed_point = false,
            jacobian_condition = NaN, smallest_singular_value = NaN,
            generalized_alpha = NaN, descriptor_alpha_difference = NaN,
            finite_generalized_count = 0, reference_alpha = NaN, alpha_difference = NaN)
    try
        nw = build_confirmatory_network(parse_portfolio(portfolio); controller_delta = c, bounds = :discovery)
        eq = initialize_with_tolerance(nw, 1e-10)
        sys = linearize_network(eq.state)
        red = reduce_dae(sys)
        cond = PD39._conditioning(sys, red)
        ref = reduced_spectrum(red)
        common = (status = "ok", equilibrium_residual = eq.residual, fixed_point = eq.fixed_point,
                  jacobian_condition = cond.jacobian_condition,
                  smallest_singular_value = cond.smallest_singular_value,
                  reference_alpha = ref.alpha)

        values = Tuple{String,String,Any}[("eigen", "Float64", eigen(red.A).values),
                                          ("eigvals", "Float64", eigvals(red.A)),
                                          ("complex_eigen", "ComplexF64", eigen(ComplexF64.(red.A)).values)]
        if precision
            push!(values, ("eigen", "Float32", eigen(Float32.(red.A)).values))
            push!(values, ("eigen", "BigFloat_matrix", eigen(BigFloat.(red.A)).values))
        end
        for (method, arithmetic, vals) in values
            λ = ComplexF64.(vals)
            keep = findall(abs.(λ) .> 1e-8)
            k = keep[argmax(real.(λ[keep]))]
            push!(rows, merge(base, common, (status = "ok", eig_method = method,
                arithmetic = arithmetic, alpha = real(λ[k]), critical_real = real(λ[k]),
                critical_imag = imag(λ[k]), alpha_difference = real(λ[k]) - ref.alpha)))
        end

        descriptor_status = "ok"
        generalized_alpha = NaN
        difference = NaN
        finite_count = 0
        descriptor_error = ""
        try
            gm = eigvals(Matrix(sys.A), Matrix(sys.M))
            finite = ComplexF64[x for x in gm if isfinite(real(x)) && isfinite(imag(x)) && abs(x) > 1e-8]
            finite_count = length(finite)
            generalized_alpha = maximum(real.(finite))
            difference = generalized_alpha - ref.alpha
        catch err
            descriptor_status = "unavailable"
            descriptor_error = sprint(showerror, err)
        end
        push!(rows, merge(base, common, (audit_type = "descriptor_crosscheck",
            status = descriptor_status, error_type = descriptor_status == "ok" ? "" : "api_or_matrix_failure",
            error_message = descriptor_error, eig_method = "generalized_eigvals_A_M",
            arithmetic = "Float64", alpha = generalized_alpha, critical_real = generalized_alpha,
            generalized_alpha = generalized_alpha, descriptor_alpha_difference = difference,
            finite_generalized_count = finite_count, alpha_difference = difference)))
    catch err
        push!(rows, merge(base, (error_type = "exception", error_message = sprint(showerror, err))))
    end
    return rows
end

function finite_perturbation_rows()
    perturbations = [
        ("pll", :control, 1, 0.005), ("filter", :control, 2, 0.005),
        ("current_control", :control, 3, 0.005), ("load_P", :load, 1, 0.005),
        ("load_Q", :load, 2, 0.005), ("IBR_P", :ibr, 1, 0.005),
        ("line_1", :branch, 1, 0.01), ("line_46", :branch, 46, 0.01),
    ]
    rows = NamedTuple[]
    for scenario in (SCENARIOS[1], only(filter(s -> s.id == HIGH_PLL_SCENARIO, SCENARIOS)))
        c0 = scenario_controller(scenario)
        ref = try
            run_margin_case(build_confirmatory_network(CANDIDATE_SG_BUSES;
                controller_delta = c0, bounds = :discovery))
        catch err
            (status = "failed", alpha = NaN, margin = NaN, stable = false,
             error_type = "reference_exception", error_message = sprint(showerror, err),
             equilibrium_residual = NaN)
        end
        for (label, kind, idx, magnitude) in perturbations, sign in (-1.0, 1.0)
            cd = c0
            ld = (0.0, 0.0)
            idelta = 0.0
            bd = zeros(Float64, 46)
            actual = sign * magnitude
            if kind == :control
                cc = collect(c0); cc[idx] += actual; cd = Tuple(cc)
            elseif kind == :load
                ll = collect(ld); ll[idx] = actual; ld = Tuple(ll)
            elseif kind == :ibr
                idelta = actual
            elseif kind == :branch
                bd[idx] = actual
            end
            result = try
                run_margin_case(build_confirmatory_network(CANDIDATE_SG_BUSES;
                    controller_delta = cd, load_delta = ld, ibr_delta = idelta,
                    branch_delta = bd, bounds = :discovery))
            catch err
                (status = "failed", alpha = NaN, margin = NaN, stable = false,
                 error_type = "exception", error_message = sprint(showerror, err),
                 equilibrium_residual = NaN)
            end
            push!(rows, (portfolio = V8, scenario = scenario.id, coordinate = label,
                sign = Int(sign), normalized_delta = actual, status = String(result.status),
                error_type = String(result.error_type), error_message = String(result.error_message),
                alpha = Float64(result.alpha), margin = Float64(result.margin),
                stable = Bool(result.stable), equilibrium_residual = Float64(result.equilibrium_residual),
                reference_alpha = Float64(ref.alpha), delta_alpha = Float64(result.alpha) - Float64(ref.alpha),
                reference_status = String(ref.status)))
        end
    end
    return DataFrame(rows)
end

selection = selection_rows()
CSV.write(joinpath(RESULTS, "PD39_COMPOSITION_SELECTION.csv"), selection)

composition = reuse_selected_composition(selection)
CSV.write(joinpath(RESULTS, "PD39_PENETRATION_COMPOSITION.csv"), composition)

pairwise = isempty(FROZEN_PAIRS) ? DataFrame() : vcat([pair_rows(p, composition) for p in FROZEN_PAIRS]...)
CSV.write(joinpath(RESULTS, "PD39_PENETRATION_COMPOSITION_PAIRS.csv"), pairwise)
CSV.write(joinpath(RESULTS, "PD39_PENETRATION_COMPOSITION_PAIR_SUMMARY.csv"), pair_summaries(pairwise))

numeric_portfolios = unique(vcat([V8], [FROZEN_PAIRS[8].a, FROZEN_PAIRS[8].b,
                                         FROZEN_PAIRS[10].a, FROZEN_PAIRS[10].b]))
numeric_cases = [s for s in SCENARIOS if s.id in ("nominal", HIGH_PLL_SCENARIO)]

tolerance_rows = NamedTuple[]
for p in [V8, FROZEN_PAIRS[8].a, FROZEN_PAIRS[8].b]
    for s in numeric_cases, tol in TOLERANCES
        println("tolerance ", p, " / ", s.id, " / ", tol)
        flush(stdout)
        push!(tolerance_rows, tolerance_audit_row(p, s, tol))
    end
end
CSV.write(joinpath(RESULTS, "PD39_NUMERICAL_TOLERANCE_SWEEP.csv"), DataFrame(tolerance_rows))

spectrum = NamedTuple[]
for p in numeric_portfolios, s in numeric_cases
    println("spectrum ", p, " / ", s.id)
    flush(stdout)
    precision = p == V8 || p == FROZEN_PAIRS[8].a || p == FROZEN_PAIRS[8].b
    append!(spectrum, spectrum_rows(p, s; precision))
end
CSV.write(joinpath(RESULTS, "PD39_NUMERICAL_EIGENSOLVER_AUDIT.csv"), DataFrame(spectrum))

perturbations = finite_perturbation_rows()
CSV.write(joinpath(RESULTS, "PD39_V8_FINITE_PERTURBATIONS.csv"), perturbations)

println("follow-up audit complete")
println("selection rows = ", nrow(selection))
println("composition rows = ", nrow(composition))
println("pairwise rows = ", nrow(pairwise))
println("tolerance rows = ", length(tolerance_rows))
println("eigensolver rows = ", length(spectrum))
println("finite perturbation rows = ", nrow(perturbations))
