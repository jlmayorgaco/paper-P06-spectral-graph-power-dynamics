using CSV
using DataFrames
using Dates
using LinearAlgebra
using NetworkDynamics
using PowerDynamics
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const RESULTS = joinpath(ROOT, "results")
mkpath(RESULTS)
const DISCOVERY = CSV.read(joinpath(RESULTS, "pd39", "portfolio_campaign", "portfolio_scenario_results.csv"), DataFrame)
const STATIC = CSV.read(joinpath(RESULTS, "PD39_CLASSICAL_BASELINES.csv"), DataFrame)
const SCENARIOS = uncertainty_scenarios()
const HIGH_PLL = filter(s -> s.pll_scale == 1.2, SCENARIOS)
const V8 = portfolio_string(CANDIDATE_SG_BUSES)
const PREDECESSORS = Dict(portfolio_string(filter(!=(b), CANDIDATE_SG_BUSES)) => b for b in CANDIDATE_SG_BUSES)
const AUDIT_PORTFOLIOS = vcat([V8], sort(collect(keys(PREDECESSORS))))
const TOLERANCES = [1e-8, 1e-10, 1e-12]
const HIGH_PLL_SCENARIO = "pll1.2_xf0.8_cc1.2"

parse_portfolio(s::AbstractString) = isempty(s) || s == "none" ? Int[] : parse.(Int, split(s, ";"))
float_or_nan(x) = x isa Missing ? NaN : Float64(x)
text_or_empty(x) = x isa Missing ? "" : String(x)
static_row(p) = only(filter(r -> String(r.portfolio) == p, eachrow(STATIC)))

function scenario_controller(s)
    return (s.pll_scale - 1, s.filter_scale - 1, s.current_control_scale - 1)
end

function frozen_discovery_row(p, sid)
    rows = filter(r -> String(r.portfolio) == p && String(r.scenario) == sid, eachrow(DISCOVERY))
    isempty(rows) ? nothing : only(rows)
end

function eq_residual(state)
    nw = extract_nw(state)
    du = zeros(Float64, length(uflat(state)))
    nw(du, uflat(state), pflat(state), state.t)
    return maximum(abs, du)
end

function critical_from_values(values)
    λ = ComplexF64.(values)
    keep = findall(abs.(λ) .> 1e-8)
    isempty(keep) && error("no non-gauge eigenvalues")
    k = keep[argmax(real.(λ[keep]))]
    return (alpha = real(λ[k]), critical = λ[k], count = length(λ), gauge_count = length(λ) - length(keep))
end

function solve_case(p, scenario, tol)
    nw = build_confirmatory_network(parse_portfolio(p);
        controller_delta = scenario_controller(scenario), bounds = :discovery)
    pf = solve_powerflow(nw; verbose = false, sparse = false,
        tol = tol, abstol = tol, reltol = tol)
    state = initialize_from_pf!(nw; pfs = pf, verbose = false, sparsepf = false,
        check = :none, tol = tol)
    return (state = state, pf_finite = all(isfinite, uflat(pf)),
        state_finite = all(isfinite, uflat(state)), fixed_point = isfixpoint(state; tol = tol),
        residual = eq_residual(state))
end

function finite_difference_jacobian(state; relative_step = 1e-6)
    nw = extract_nw(state)
    u = Float64.(uflat(state))
    p = Float64.(pflat(state))
    t = state.t
    n = length(u)
    J = Matrix{Float64}(undef, n, n)
    for j in 1:n
        h = relative_step * (1 + abs(u[j]))
        up = copy(u); um = copy(u)
        up[j] += h; um[j] -= h
        fp = zeros(Float64, n); fm = zeros(Float64, n)
        nw(fp, up, p, t); nw(fm, um, p, t)
        J[:, j] = (fp .- fm) ./ (2h)
    end
    return J
end

function reduce_dae_independent(A, M)
    if M isa UniformScaling
        return A
    end
    md = diag(Matrix(M))
    didx = findall(==(1), md)
    aidx = findall(==(0), md)
    isempty(aidx) && return A[didx, didx]
    return A[didx, didx] - A[didx, aidx] * (A[aidx, aidx] \ A[aidx, didx])
end

function descriptor_alpha(A, M)
    vals = eigvals(Matrix(A), Matrix(M))
    finite = ComplexF64[x for x in vals if isfinite(real(x)) && isfinite(imag(x)) && abs(x) > 1e-8]
    isempty(finite) && error("no finite non-gauge descriptor eigenvalues")
    return (alpha = maximum(real.(finite)), count = length(finite))
end

function base_a_row(p, s, tol)
    common = (portfolio = p, scenario = s.id, audit_layer = "tolerance",
        tolerance = tol, method = "reference_eigen", arithmetic = "Float64",
        status = "failed", error_type = "", error_message = "", alpha = NaN,
        critical_eigenvalue = "", critical_frequency_hz = NaN, damping_ratio = NaN,
        equilibrium_residual = NaN, fixed_point = false, jacobian_condition = NaN,
        smallest_singular_value = NaN, reference_alpha = NaN, alpha_difference = NaN,
        independent_fd_alpha = NaN, fd_alpha_difference = NaN, descriptor_alpha = NaN,
        descriptor_alpha_difference = NaN, descriptor_finite_count = 0,
        fd_descriptor_alpha = NaN, fd_descriptor_difference = NaN)
    try
        eq = solve_case(p, s, tol)
        sys = linearize_network(eq.state)
        red = reduce_dae(sys)
        sp = critical_from_values(eigen(red.A).values)
        sv = svdvals(red.A)
        cond = maximum(sv) / max(minimum(sv), eps(Float64))
        desc = try descriptor_alpha(sys.A, sys.M) catch; (alpha = NaN, count = 0) end
        return merge(common, (status = "ok", alpha = sp.alpha, critical_eigenvalue = string(sp.critical),
            critical_frequency_hz = abs(imag(sp.critical)) / (2pi),
            damping_ratio = abs(real(sp.critical)) / max(abs(sp.critical), eps(Float64)),
            equilibrium_residual = eq.residual, fixed_point = eq.fixed_point,
            jacobian_condition = cond, smallest_singular_value = minimum(sv),
            reference_alpha = sp.alpha, alpha_difference = 0.0,
            descriptor_alpha = desc.alpha, descriptor_alpha_difference = desc.alpha - sp.alpha,
            descriptor_finite_count = desc.count))
    catch err
        return merge(common, (error_type = "exception", error_message = sprint(showerror, err)))
    end
end

function numerical_truth_audit()
    rows = NamedTuple[]
    modal_rows = NamedTuple[]
    shapes = Dict{Tuple{String,String},Vector{ComplexF64}}()
    tight_states = Dict{Tuple{String,String},Any}()

    for p in AUDIT_PORTFOLIOS, s in SCENARIOS, tol in TOLERANCES
        println("A tolerance ", p, " / ", s.id, " / ", tol)
        flush(stdout)
        push!(rows, base_a_row(p, s, tol))
        if tol == 1e-10
            try
                eq = solve_case(p, s, tol)
                tight_states[(p, s.id)] = eq.state
            catch
            end
        end
    end

    # Independent central-difference Jacobian, eigensolver, and descriptor check
    # at the tight reference state for every one of the 81 required cases.
    for p in AUDIT_PORTFOLIOS, s in SCENARIOS
        println("A independent ", p, " / ", s.id)
        flush(stdout)
        state = get(tight_states, (p, s.id), nothing)
        common = (portfolio = p, scenario = s.id, audit_layer = "independent",
            tolerance = 1e-10, method = "central_finite_difference", arithmetic = "Float64",
            status = "failed", error_type = "", error_message = "", alpha = NaN,
            critical_eigenvalue = "", critical_frequency_hz = NaN, damping_ratio = NaN,
            equilibrium_residual = NaN, fixed_point = false, jacobian_condition = NaN,
            smallest_singular_value = NaN, reference_alpha = NaN, alpha_difference = NaN,
            independent_fd_alpha = NaN, fd_alpha_difference = NaN, descriptor_alpha = NaN,
            descriptor_alpha_difference = NaN, descriptor_finite_count = 0,
            fd_descriptor_alpha = NaN, fd_descriptor_difference = NaN)
        if state === nothing
            push!(rows, merge(common, (error_type = "missing_tight_state", error_message = "1e-10 initialization failed")))
            continue
        end
        try
            sys = linearize_network(state)
            red = reduce_dae(sys)
            ref = critical_from_values(eigen(red.A).values)
            Afd = finite_difference_jacobian(state)
            redfd = reduce_dae_independent(Afd, sys.M)
            fd = critical_from_values(eigvals(redfd))
            desc = descriptor_alpha(sys.A, sys.M)
            fddesc = descriptor_alpha(Afd, sys.M)
            sv = svdvals(red.A)
            cond = maximum(sv) / max(minimum(sv), eps(Float64))
            push!(rows, merge(common, (status = "ok", alpha = fd.alpha,
                critical_eigenvalue = string(fd.critical),
                critical_frequency_hz = abs(imag(fd.critical)) / (2pi),
                damping_ratio = abs(real(fd.critical)) / max(abs(fd.critical), eps(Float64)),
                equilibrium_residual = eq_residual(state), fixed_point = isfixpoint(state; tol = 1e-10),
                jacobian_condition = cond, smallest_singular_value = minimum(sv),
                reference_alpha = ref.alpha, alpha_difference = fd.alpha - ref.alpha,
                independent_fd_alpha = fd.alpha, fd_alpha_difference = fd.alpha - ref.alpha,
                descriptor_alpha = desc.alpha, descriptor_alpha_difference = desc.alpha - ref.alpha,
                descriptor_finite_count = desc.count, fd_descriptor_alpha = fddesc.alpha,
                fd_descriptor_difference = fddesc.alpha - fd.alpha)))

            m = modal_report(state)
            mac = NaN
            shape = try common_mode_shape(state) catch; nothing end
            shape !== nothing && (shapes[(p, s.id)] = shape)
            if p != V8 && haskey(shapes, (V8, s.id)) && shape !== nothing
                mac = modal_assurance(shapes[(V8, s.id)], shape)
            end
            push!(modal_rows, (portfolio = p, scenario = s.id, status = "ok",
                alpha = real(m.critical_eigenvalue), critical_eigenvalue = string(m.critical_eigenvalue),
                critical_frequency_hz = m.critical_frequency_hz, damping_ratio = m.damping_ratio,
                critical_mode_family = m.critical_mode_family, critical_state_labels = m.critical_state_labels,
                critical_participation = m.critical_participation, equilibrium_residual = m.equilibrium_residual,
                jacobian_condition = m.jacobian_condition, smallest_singular_value = m.smallest_singular_value,
                g_z_condition = m.g_z_condition, eigenvector_condition = m.eigenvector_condition,
                mac_to_v8 = mac))
        catch err
            push!(rows, merge(common, (error_type = "independent_audit_exception", error_message = sprint(showerror, err))))
            push!(modal_rows, (portfolio = p, scenario = s.id, status = "failed", alpha = NaN,
                critical_eigenvalue = "", critical_frequency_hz = NaN, damping_ratio = NaN,
                critical_mode_family = "", critical_state_labels = "", critical_participation = "",
                equilibrium_residual = NaN, jacobian_condition = NaN, smallest_singular_value = NaN,
                g_z_condition = NaN, eigenvector_condition = NaN, mac_to_v8 = NaN))
        end
    end
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv"), DataFrame(rows))
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_MODE_TRACKING.csv"), DataFrame(modal_rows))
    return (rows = DataFrame(rows), modal = DataFrame(modal_rows), states = tight_states)
end

function best_worst_7of8()
    d = combine(groupby(DISCOVERY, :portfolio), :dynamic_margin => minimum => :m9)
    d = innerjoin(d, STATIC[:, [:portfolio, :cardinality]], on = :portfolio)
    d = filter(r -> r.cardinality == 7, d)
    best = first(sort(collect(eachrow(d)), by = r -> (-r.m9, r.portfolio)))
    worst = first(sort(collect(eachrow(d)), by = r -> (r.m9, r.portfolio)))
    return (best = String(best.portfolio), worst = String(worst.portfolio), best_m9 = best.m9, worst_m9 = worst.m9)
end

function write_tds_trace(label, scenario, td)
    path = joinpath(RESULTS, "PD39_255PLUS1_TDS_TRACE_$(label)_$(scenario.id).csv")
    CSV.write(path, DataFrame(time = td.t, frequency_spread_hz = td.frequency_spread_hz,
        voltage_pu = td.voltage, voltage_deviation_pu = td.voltage_deviation_pu))
    return path
end

function high_pll_tds(audit)
    bw = best_worst_7of8()
    cases = [("V8", V8), ("best_7of8", bw.best), ("worst_7of8", bw.worst)]
    rows = NamedTuple[]
    for (label, p) in cases, s in HIGH_PLL
        println("B TDS ", label, " / ", s.id)
        flush(stdout)
        alpha_rows = filter(r -> r.portfolio == p && r.scenario == s.id && r.audit_layer == "independent",
                            eachrow(audit.rows))
        linear_alpha = isempty(alpha_rows) ? NaN : Float64(first(alpha_rows).independent_fd_alpha)
        try
            td = simulate_tds_case(parse_portfolio(p); pulse = 0.01,
                controller_delta = scenario_controller(s), bounds = :discovery,
                tspan = (0.0, 20.0), saveat = 0.01)
            trace = write_tds_trace(label, s, td)
            sign_consistent = isfinite(linear_alpha) && isfinite(td.estimated_rate_s_inv) &&
                sign(linear_alpha) == sign(td.estimated_rate_s_inv)
            push!(rows, (case_id = label, portfolio = p, scenario = s.id, pulse = 0.01,
                status = "ok", error = "", linear_alpha = linear_alpha,
                estimated_rate_s_inv = td.estimated_rate_s_inv,
                sign_consistent = sign_consistent,
                max_frequency_spread_hz = td.max_frequency_spread_hz,
                max_voltage_deviation_pu = td.max_voltage_deviation_pu,
                frequency_settling_s = td.frequency_settling_s,
                voltage_settling_s = td.voltage_settling_s,
                voltage_min_pu = td.voltage_min_pu, voltage_max_pu = td.voltage_max_pu,
                trace = trace))
        catch err
            push!(rows, (case_id = label, portfolio = p, scenario = s.id, pulse = 0.01,
                status = "failed", error = sprint(showerror, err), linear_alpha = linear_alpha,
                estimated_rate_s_inv = NaN, sign_consistent = false,
                max_frequency_spread_hz = NaN, max_voltage_deviation_pu = NaN,
                frequency_settling_s = NaN, voltage_settling_s = NaN,
                voltage_min_pu = NaN, voltage_max_pu = NaN, trace = ""))
        end
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_HIGH_PLL_TDS.csv"), out)
    return (table = out, best_worst = bw)
end

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

function composition_mechanism(audit)
    rows = NamedTuple[]
    for pair in COMPOSITION_PAIRS, s in SCENARIOS
        ra = filter(r -> r.portfolio == pair.a && r.scenario == s.id && r.audit_layer == "independent",
                    eachrow(audit.rows))
        rb = filter(r -> r.portfolio == pair.b && r.scenario == s.id && r.audit_layer == "independent",
                    eachrow(audit.rows))
        sa, sb = static_row(pair.a), static_row(pair.b)
        ok = !isempty(ra) && !isempty(rb) && first(ra).status == "ok" && first(rb).status == "ok"
        push!(rows, (pair_id = pair.pair_id, cardinality = pair.cardinality, scenario = s.id,
            portfolio_a = pair.a, portfolio_b = pair.b,
            delta_converted_mw = abs(Float64(sa.converted_mw) - Float64(sb.converted_mw)),
            delta_ibr_mva = abs(Float64(sa.ibr_mva) - Float64(sb.ibr_mva)),
            delta_remaining_sg_mw = abs(Float64(sa.remaining_sg_mw) - Float64(sb.remaining_sg_mw)),
            delta_remaining_sg_mva = abs(Float64(sa.remaining_sg_mva) - Float64(sb.remaining_sg_mva)),
            delta_remaining_inertia_mva_s = abs(Float64(sa.remaining_inertia_mva_s) - Float64(sb.remaining_inertia_mva_s)),
            alpha_a = ok ? first(ra).independent_fd_alpha : NaN,
            alpha_b = ok ? first(rb).independent_fd_alpha : NaN,
            delta_alpha_a_minus_b = ok ? first(ra).independent_fd_alpha - first(rb).independent_fd_alpha : NaN,
            mode_mac_a_to_v8 = ok ? first(ra).mac_to_v8 : NaN,
            mode_mac_b_to_v8 = ok ? first(rb).mac_to_v8 : NaN,
            status = ok ? "ok" : "missing_mode_audit"))
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_PENETRATION_COMPOSITION.csv"), out)

    source_model = read(joinpath(ROOT, "src", "pd39", "model.jl"), String)
    source_confirmatory = read(joinpath(ROOT, "src", "pd39", "confirmatory.jl"), String)
    has_continuous = occursin(r"(?i)homotopy|interpolat|blend|morph", source_model * source_confirmatory)
    homotopy = DataFrame((item = ["SG_to_GFL_physical_homotopy"],
        status = [has_continuous ? "requires_manual_semantic_review" : "BLOCKED"],
        documented_continuous_path = [has_continuous],
        reason = [has_continuous ? "candidate interpolation token found; no automatic homotopy run" :
            "model exposes discrete replace_buses/compile_bus replacement and no documented continuous SG-GFL interpolation"],
        action = ["No artificial blend constructed"]))
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_HOMOTOPY_AUDIT.csv"), homotopy)
    return out, homotopy
end

function closure_audit()
    relevant = String[]
    for (dir, _, files) in walkdir(ROOT)
        occursin(".git", dir) && continue
        for f in files
            (endswith(f, ".jl") || endswith(f, ".md") || endswith(f, ".toml")) || continue
            path = joinpath(dir, f)
            text = lowercase(read(path, String))
            if occursin("determinant identity", text) || occursin("network closure", text) ||
               (occursin("q→-1", text) || occursin("q->-1", text))
                push!(relevant, path)
            end
        end
    end
    out = DataFrame(item = ["K", "D", "Q", "determinant_identity", "Q_to_minus_one_boundary"],
        status = fill("BLOCKED", 5),
        constructible = fill(false, 5),
        evidence = fill("No documented exact K,D,Q closure objects, dimensions, units, or identity in the installed PD39 model/API; no proxy constructed.", 5),
        repository_hits = fill(join(relevant, "|"), 5))
    CSV.write(joinpath(RESULTS, "PD39_255PLUS1_CLOSURE_AUDIT.csv"), out)
    return out
end

println("PD39 255+1 mechanism validation start")
audit = numerical_truth_audit()
tds = high_pll_tds(audit)
composition, homotopy = composition_mechanism(audit)
closure = closure_audit()
println("A rows=", nrow(audit.rows), " modal=", nrow(audit.modal))
println("B TDS rows=", nrow(tds.table))
println("D composition rows=", nrow(composition), " homotopy=", homotopy.status[1])
println("E closure=", unique(closure.status))
println("PD39 255+1 mechanism validation complete")
