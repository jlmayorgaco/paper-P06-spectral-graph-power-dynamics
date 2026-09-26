using CSV
using DataFrames
using Graphs
using LinearAlgebra
using NetworkDynamics
using PowerDynamics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS = joinpath(ROOT, "results")
mkpath(RESULTS)
const CENSUS = CSV.read(joinpath(RESULTS, "PD39_255PLUS1_HOLDOUT_CENSUS.csv"), DataFrame)
const STATIC = CSV.read(joinpath(RESULTS, "PD39_CLASSICAL_BASELINES.csv"), DataFrame)
const CONDITIONS = CSV.read(joinpath(RESULTS, "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
const AUDIT_CONDITIONS = ["C12", "nominal"]
const TOLERANCES = [1e-8, 1e-10, 1e-12]
const FD_STEPS = [1e-5, 1e-6, 1e-7]
const SELECTED = [
    (blocker_id = "H1", condition_id = "C12", portfolio = "35;36", control = "32;36"),
    (blocker_id = "H2", condition_id = "C12", portfolio = "37;38", control = "36;38"),
    (blocker_id = "H3", condition_id = "C12", portfolio = "30;32;33", control = "30;33;35"),
]

parse_portfolio(p::AbstractString) = p == "none" ? Int[] : parse.(Int, split(p, ";"))
portfolio_string(p) = isempty(p) ? "none" : join(sort(Int.(collect(p))), ";")
static_row(p) = only(filter(r -> String(r.portfolio) == p, eachrow(STATIC)))
census_row(p, c) = begin
    rs = filter(r -> String(r.portfolio) == p && String(r.condition) == c, eachrow(CENSUS))
    isempty(rs) ? nothing : only(rs)
end

function condition_row(id)
    id == "nominal" && return (id = id, delta_pll = 0.0, delta_xf = 0.0, delta_cc = 0.0,
        delta_load_p = 0.0, delta_load_q = 0.0, delta_ibr_p = 0.0, branches = zeros(Float64, 46))
    r = only(filter(x -> String(x.condition) == id, eachrow(CONDITIONS)))
    b = [Float64(r[Symbol("delta_line_$(lpad(i, 2, '0'))")]) for i in 1:46]
    return (id = id, delta_pll = Float64(r.delta_pll), delta_xf = Float64(r.delta_xf),
        delta_cc = Float64(r.delta_cc), delta_load_p = Float64(r.delta_load_p),
        delta_load_q = Float64(r.delta_load_q), delta_ibr_p = Float64(r.delta_ibr_p), branches = b)
end

scenario_kwargs(c) = (controller_delta = (c.delta_pll, c.delta_xf, c.delta_cc),
    load_delta = (c.delta_load_p, c.delta_load_q), ibr_delta = c.delta_ibr_p,
    branch_delta = c.branches)

function proper_subsets(p)
    buses = parse_portfolio(p)
    out = String[]
    for mask in 0:(2^length(buses) - 2)
        push!(out, portfolio_string([buses[i] for i in eachindex(buses) if (mask >> (i - 1)) & 1 == 1]))
    end
    unique(out)
end

function selected_roles()
    roles = Dict{String,Vector{String}}()
    function add!(p, x)
        roles[p] = unique(vcat(get(roles, p, String[]), [x]))
    end
    for h in SELECTED
        add!(h.portfolio, "selected_blocker_$(h.blocker_id)")
        for s in proper_subsets(h.portfolio)
            add!(s, "proper_subset_of_$(h.blocker_id)")
        end
        add!(h.control, "stable_matched_control_for_$(h.blocker_id)")
    end
    roles
end

function write_selection_skeleton()
    rows = NamedTuple[]
    for h in SELECTED
        bm = static_row(h.portfolio); cm = static_row(h.control)
        br = census_row(h.portfolio, h.condition_id); cr = census_row(h.control, h.condition_id)
        push!(rows, (blocker_id = h.blocker_id, condition_id = h.condition_id, portfolio = h.portfolio,
            cardinality = Int(bm.cardinality), converted_MW = Float64(bm.converted_mw),
            converted_MVA = Float64(bm.ibr_mva), remaining_inertia = Float64(bm.remaining_inertia_mva_s),
            alpha = br === nothing ? NaN : Float64(br.alpha), frequency = NaN,
            stable_control = h.control, control_MW = Float64(cm.converted_mw),
            control_MVA = Float64(cm.ibr_mva), control_inertia = Float64(cm.remaining_inertia_mva_s),
            control_alpha = cr === nothing ? NaN : Float64(cr.alpha), control_frequency = NaN,
            aggregate_distance = NaN, selection_source = "frozen_C12_census_before_new_numerics"))
    end
    CSV.write(joinpath(RESULTS, "PD39_SELECTED_BLOCKERS.csv"), DataFrame(rows))
end

function eq_residual(state)
    nw = extract_nw(state)
    du = zeros(Float64, length(uflat(state)))
    nw(du, uflat(state), pflat(state), state.t)
    maximum(abs, du)
end

function voltage_range(state)
    vals = Float64[]
    for bus in 1:39
        try
            ur = Float64(state[VIndex(bus, :busbar₊u_r)])
            ui = Float64(state[VIndex(bus, :busbar₊u_i)])
            isfinite(ur) && isfinite(ui) && push!(vals, hypot(ur, ui))
        catch
        end
    end
    isempty(vals) ? (minimum = NaN, maximum = NaN) : (minimum = minimum(vals), maximum = maximum(vals))
end

function solve_case(p, c, tol)
    nw = build_confirmatory_network(parse_portfolio(p); scenario_kwargs(c)..., bounds = :primary)
    pf = solve_powerflow(nw; verbose = false, sparse = false, tol = tol, abstol = tol, reltol = tol)
    state = initialize_from_pf!(nw; pfs = pf, verbose = false, sparsepf = false, check = :none, tol = tol)
    (nw = nw, pf = pf, state = state, pf_finite = all(isfinite, uflat(pf)),
        state_finite = all(isfinite, uflat(state)), fixed_point = isfixpoint(state; tol = tol),
        residual = eq_residual(state))
end

function critical(values; gauge_tol = 1e-8)
    λ = ComplexF64.(values)
    finite = filter(x -> isfinite(real(x)) && isfinite(imag(x)), λ)
    keep = filter(x -> abs(x) > gauge_tol, finite)
    isempty(keep) && error("no finite non-gauge eigenvalue")
    z = keep[argmax(real.(keep))]
    (alpha = real(z), value = z, frequency = abs(imag(z)) / (2pi),
        damping = abs(real(z)) / max(abs(z), eps(Float64)), finite_count = length(finite),
        gauge_count = length(finite) - length(keep), nearest_zero = minimum(abs.(keep)))
end

function independent_reduce(A, M)
    M isa UniformScaling && return Matrix(A)
    md = diag(Matrix(M))
    didx = findall(x -> abs(x) > 0.5, md); aidx = findall(x -> abs(x) <= 0.5, md)
    isempty(aidx) && return Matrix(A)[didx, didx]
    AA = Matrix(A)
    AA[didx, didx] - AA[didx, aidx] * (AA[aidx, aidx] \ AA[aidx, didx])
end

function descriptor_critical(A, M)
    M isa UniformScaling && return (critical(Matrix(A))..., status = "uniform_scaling_equivalent")
    (critical(eigvals(Matrix(A), Matrix(M)))..., status = "generalized_dense")
end

function fd_jacobian(state, relative_step)
    nw = extract_nw(state); u = Float64.(uflat(state)); p = Float64.(pflat(state)); t = state.t
    n = length(u); J = Matrix{Float64}(undef, n, n)
    for j in 1:n
        h = relative_step * (1 + abs(u[j]))
        up = copy(u); um = copy(u); up[j] += h; um[j] -= h
        fp = zeros(Float64, n); fm = zeros(Float64, n)
        nw(fp, up, p, t); nw(fm, um, p, t)
        J[:, j] = (fp .- fm) ./ (2h)
    end
    J
end

function classify_mode(freq, participation, gz_condition)
    !isfinite(gz_condition) || gz_condition > 1e12 && return "DAE_SINGULARITY"
    entries = split(String(participation), "|"); sync = 0.0; conv = 0.0
    for item in entries
        parts = split(item, ":"); length(parts) < 2 && continue
        score = try parse(Float64, parts[end]) catch; 0.0 end
        label = join(parts[1:end-1], ":")
        occursin(r"machine|gov|avr|ω|delta|δ", label) && (sync += score)
        occursin(r"SimpleGFLDC|PLL|CC1|gfl|dc", label) && (conv += score)
    end
    if freq >= 0.1 && conv >= 0.2 && sync >= 0.2
        "MIXED_OSCILLATORY"
    elseif freq >= 0.1 && conv >= 0.25
        "CONVERTER_OSCILLATORY"
    elseif freq >= 0.1 && sync >= 0.25
        "EM_OSCILLATORY"
    elseif freq < 0.1 && (conv + sync) >= 0.25
        "APERIODIC_DYNAMIC"
    else
        "NUMERICALLY_UNRESOLVED"
    end
end

function audit_one(p, c, tol; do_fd = false)
    base = (portfolio = p, condition = c.id, tolerance = tol, status = "failed", error_type = "",
        error_message = "", equilibrium_converged = false, equilibrium_residual = NaN,
        voltage_min_pu = NaN, voltage_max_pu = NaN, load_P_MW = NaN, load_Q_MVAr = NaN,
        requested_ibr_P_MW = NaN, actual_injection_status = "not_exposed_by_stock_api",
        actual_P_MW = NaN, actual_Q_MVAr = NaN, power_balance_residual = NaN, load_delivered = false,
        native_alpha = NaN, independent_alpha = NaN, descriptor_alpha = NaN,
        native_frequency_hz = NaN, independent_frequency_hz = NaN, descriptor_frequency_hz = NaN,
        native_critical_eigenvalue = "", independent_critical_eigenvalue = "",
        descriptor_critical_eigenvalue = "", native_dense_alpha_difference = NaN,
        descriptor_alpha_difference = NaN, jacobian_condition = NaN, smallest_singular_value = NaN,
        g_z_condition = NaN, g_z_smallest_singular_value = NaN, eigenvector_condition = NaN,
        nearest_nonzero_abs = NaN, gauge_count = 0, critical_mode_family = "",
        mechanism_class = "", critical_participation = "", fd_alpha_1e5 = NaN,
        fd_alpha_1e6 = NaN, fd_alpha_1e7 = NaN, fd_difference_1e5 = NaN,
        fd_difference_1e6 = NaN, fd_difference_1e7 = NaN, independent_sign_agreement = false,
        descriptor_sign_agreement = false, fd_sign_agreement = false)
    try
        eq = solve_case(p, c, tol); vr = voltage_range(eq.state); data = ieee39_data()
        loadp = sum(Float64(r.Pset) for r in eachrow(data.load)) * (1 + c.delta_load_p) * 100
        loadq = sum(Float64(r.Qset) for r in eachrow(data.load)) * (1 + c.delta_load_q) * 100
        ibrp = sum(Float64(data.bus[data.bus.bus .== b, :P][1]) for b in parse_portfolio(p)) *
            (1 + c.delta_ibr_p) * 100
        if !eq.pf_finite || !eq.state_finite || !eq.fixed_point
            return merge(base, (error_type = "equilibrium_not_qualified",
                error_message = "PF/dynamic state/fixed-point gate failed",
                equilibrium_residual = eq.residual, voltage_min_pu = vr.minimum,
                voltage_max_pu = vr.maximum, load_P_MW = loadp, load_Q_MVAr = loadq,
                requested_ibr_P_MW = ibrp, power_balance_residual = eq.residual))
        end
        sys = linearize_network(eq.state); red = reduce_dae(sys)
        native = stability_audit(eq.state; gauge_tol = 1e-8); dense = critical(eigvals(Matrix(red.A)))
        desc = descriptor_critical(sys.A, sys.M); sv = svdvals(Matrix(red.A)); m = modal_report(eq.state)
        fdvals = Dict{Float64,Float64}()
        if do_fd
            for h in FD_STEPS
                fdred = independent_reduce(fd_jacobian(eq.state, h), sys.M)
                fdvals[h] = critical(eigvals(fdred)).alpha
            end
        end
        return merge(base, (status = "ok", equilibrium_converged = true, equilibrium_residual = eq.residual,
            voltage_min_pu = vr.minimum, voltage_max_pu = vr.maximum, load_P_MW = loadp,
            load_Q_MVAr = loadq, requested_ibr_P_MW = ibrp, power_balance_residual = eq.residual,
            load_delivered = isfinite(vr.minimum) && vr.minimum >= 0.90 && vr.maximum <= 1.10,
            native_alpha = native.max_real, independent_alpha = dense.alpha, descriptor_alpha = desc.alpha,
            native_frequency_hz = m.critical_frequency_hz, independent_frequency_hz = dense.frequency,
            descriptor_frequency_hz = desc.frequency, native_critical_eigenvalue = string(m.critical_eigenvalue),
            independent_critical_eigenvalue = string(dense.value), descriptor_critical_eigenvalue = string(desc.value),
            native_dense_alpha_difference = dense.alpha - native.max_real,
            descriptor_alpha_difference = desc.alpha - native.max_real, jacobian_condition = m.jacobian_condition,
            smallest_singular_value = m.smallest_singular_value, g_z_condition = m.g_z_condition,
            g_z_smallest_singular_value = m.g_z_smallest_singular_value, eigenvector_condition = m.eigenvector_condition,
            nearest_nonzero_abs = dense.nearest_zero, gauge_count = dense.gauge_count,
            critical_mode_family = m.critical_mode_family,
            mechanism_class = classify_mode(m.critical_frequency_hz, m.critical_participation, m.g_z_condition),
            critical_participation = m.critical_participation,
            fd_alpha_1e5 = get(fdvals, 1e-5, NaN), fd_alpha_1e6 = get(fdvals, 1e-6, NaN),
            fd_alpha_1e7 = get(fdvals, 1e-7, NaN), fd_difference_1e5 = get(fdvals, 1e-5, NaN) - native.max_real,
            fd_difference_1e6 = get(fdvals, 1e-6, NaN) - native.max_real,
            fd_difference_1e7 = get(fdvals, 1e-7, NaN) - native.max_real,
            independent_sign_agreement = sign(native.max_real) == sign(dense.alpha),
            descriptor_sign_agreement = sign(native.max_real) == sign(desc.alpha),
            fd_sign_agreement = isfinite(get(fdvals, 1e-6, NaN)) &&
                sign(native.max_real) == sign(get(fdvals, 1e-6, NaN))))
    catch err
        merge(base, (error_type = "audit_exception", error_message = sprint(showerror, err)))
    end
end

function write_selection_final!(audit)
    rows = NamedTuple[]
    for h in SELECTED
        br = filter(r -> r.portfolio == h.portfolio && r.condition == h.condition_id && r.tolerance == 1e-10, eachrow(audit))
        cr = filter(r -> r.portfolio == h.control && r.condition == h.condition_id && r.tolerance == 1e-10, eachrow(audit))
        bm = static_row(h.portfolio); cm = static_row(h.control)
        b = isempty(br) ? nothing : only(br); cc = isempty(cr) ? nothing : only(cr)
        distance = cc === nothing ? NaN : sqrt(
            ((Float64(cm.converted_mw)-Float64(bm.converted_mw))/max(Float64(bm.converted_mw),1.0))^2 +
            ((Float64(cm.ibr_mva)-Float64(bm.ibr_mva))/max(Float64(bm.ibr_mva),1.0))^2 +
            ((Float64(cm.remaining_inertia_mva_s)-Float64(bm.remaining_inertia_mva_s))/
                max(Float64(bm.remaining_inertia_mva_s),1.0))^2)
        push!(rows, (blocker_id = h.blocker_id, condition_id = h.condition_id, portfolio = h.portfolio,
            cardinality = Int(bm.cardinality), converted_MW = Float64(bm.converted_mw),
            converted_MVA = Float64(bm.ibr_mva), remaining_inertia = Float64(bm.remaining_inertia_mva_s),
            alpha = b === nothing ? NaN : b.native_alpha, frequency = b === nothing ? NaN : b.native_frequency_hz,
            stable_control = h.control, control_MW = Float64(cm.converted_mw),
            control_MVA = Float64(cm.ibr_mva), control_inertia = Float64(cm.remaining_inertia_mva_s),
            control_alpha = cc === nothing ? NaN : cc.native_alpha,
            control_frequency = cc === nothing ? NaN : cc.native_frequency_hz,
            aggregate_distance = distance, selection_source = "frozen_selection_GateA_1e-10"))
    end
    CSV.write(joinpath(RESULTS, "PD39_SELECTED_BLOCKERS.csv"), DataFrame(rows))
end

function run_gate_a()
    write_selection_skeleton(); roles = selected_roles(); rows = NamedTuple[]
    for p in sort(collect(keys(roles))), cid in AUDIT_CONDITIONS, tol in TOLERANCES
        do_fd = tol == 1e-10; println("A ", p, " / ", cid, " / ", tol, " fd=", do_fd); flush(stdout)
        push!(rows, merge(audit_one(p, condition_row(cid), tol; do_fd = do_fd),
            (roles = join(roles[p], "|"), audit_layer = "tolerance")))
        CSV.write(joinpath(RESULTS, "PD39_BLOCKER_NUMERICAL_AUDIT.csv"), DataFrame(rows))
    end
    out = DataFrame(rows); write_selection_final!(out); out
end

println("PD39 blocker atlas Gate A start")
a = run_gate_a()
println("Gate A rows=", nrow(a), " portfolios=", length(unique(a.portfolio)),
    " conditions=", length(unique(a.condition)))
println("PD39 blocker atlas Gate A complete")
