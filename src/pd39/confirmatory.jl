using CSV
using DataFrames
using Graphs
using LinearAlgebra
using NetworkDynamics
using OrdinaryDiffEqRosenbrock
using PowerDynamics
using Random
using SciMLBase
using Statistics

export ConfirmatoryCase, build_confirmatory_network, run_confirmatory_case,
    modal_report, diagnostic_report, maximin_lhs, scenario_vector,
    candidate_dispatch_mw, portfolio_string, apply_line_gamma,
    common_mode_shape, modal_assurance, run_margin_case,
    build_tds_network, simulate_tds_case, largest_load_bus

const PRIMARY_CONTROL_BOUNDS = 0.10
const PRIMARY_LOAD_BOUNDS = 0.05
const PRIMARY_IBR_BOUNDS = 0.05
const PRIMARY_NETWORK_BOUNDS = 0.10
const CONFIRMATORY_SEED = 39024
const TDS_LOAD_BUS = 39
const VOLTAGE_MIN_PU = 0.90
const VOLTAGE_MAX_PU = 1.10

portfolio_string(buses) = isempty(buses) ? "none" : join(sort(Int.(collect(buses))), ";")

function candidate_dispatch_mw()
    d = ieee39_data()
    return DataFrame(
        bus = CANDIDATE_SG_BUSES,
        dispatch_mw = [100.0 * Float64(d.bus[findfirst(d.bus.bus .== b), :P]) for b in CANDIDATE_SG_BUSES],
    )
end

"Return an actual perturbation vector with the frozen seven global coordinates." 
function scenario_vector(control = (0.0, 0.0, 0.0), operating = (0.0, 0.0, 0.0),
                         branches = zeros(Float64, 46))
    return vcat(Float64.(collect(control)), Float64.(collect(operating)), Float64.(branches))
end

"Construct a maximin Latin-hypercube in [0,1]^d with a deterministic seed." 
function maximin_lhs(n::Integer, d::Integer; seed::Integer = CONFIRMATORY_SEED,
                     candidates::Integer = 128)
    n > 0 && d > 0 || throw(ArgumentError("n and d must be positive"))
    rng = MersenneTwister(seed)
    best = nothing
    best_score = -Inf
    for _ in 1:candidates
        x = Matrix{Float64}(undef, n, d)
        for j in 1:d
            x[:, j] = (randperm(rng, n) .- rand(rng, n)) ./ n
        end
        score = Inf
        for i in 1:n-1, k in i+1:n
            score = min(score, sum(abs2, @view(x[i, :]) .- @view(x[k, :])))
        end
        if score > best_score
            best_score = score
            best = x
        end
    end
    return best
end

"Apply a physically normalized branch impedance perturbation to copied edges." 
function apply_line_gamma(edges, branch_delta)
    length(branch_delta) == length(edges) || throw(ArgumentError("one branch coordinate per edge required"))
    for (i, edge) in enumerate(edges)
        δ = Float64(branch_delta[i])
        isfinite(δ) || throw(ArgumentError("nonfinite branch coordinate"))
        r0 = get_default(edge, :piline₊R)
        x0 = get_default(edge, :piline₊X)
        set_default!(edge, :piline₊R, r0 * (1 + δ))
        set_default!(edge, :piline₊X, x0 * (1 + δ))
    end
    return edges
end

function _set_load_operating_point!(vertices, load_p_delta, load_q_delta)
    data = ieee39_data()
    for row in eachrow(data.load)
        bus = Int(row.bus)
        old_p = Float64(row.Pset)
        old_q = Float64(row.Qset)
        new_p = old_p * (1 + load_p_delta)
        new_q = old_q * (1 + load_q_delta)
        set_default!(vertices[bus], Regex(raw"Pset$"), new_p)
        set_default!(vertices[bus], Regex(raw"Qset$"), new_q)
        busrow = data.bus[findfirst(data.bus.bus .== bus), :]
        if String(busrow.bus_type) == "PQ"
            # For a pure load bus, the PF net injection is the ZIP operating point.
            set_pfmodel!(vertices[bus], pfPQ(P = Float64(busrow.P) * (1 + load_p_delta),
                                             Q = Float64(busrow.Q) * (1 + load_q_delta)))
        elseif String(busrow.bus_type) == "PV"
            # For a generator-plus-load bus, change net P by the load increment
            # while retaining its voltage setpoint and generator dispatch convention.
            p_net = Float64(busrow.P) + (new_p - old_p)
            set_pfmodel!(vertices[bus], pfPV(P = p_net, V = Float64(busrow.V)))
        end
    end
    return vertices
end

function _set_replaced_dispatch!(vertices, portfolio, ibr_delta)
    abs(ibr_delta) <= 1 || throw(ArgumentError("IBR dispatch multiplier outside physical domain"))
    data = ieee39_data()
    for bus in portfolio
        row = data.bus[findfirst(data.bus.bus .== bus), :]
        set_pfmodel!(vertices[bus], pfPV(P = Float64(row.P) * (1 + ibr_delta), V = Float64(row.V)))
    end
    return vertices
end

"Build a case using only model parameters or operating-point objects present in the stock model." 
function build_confirmatory_network(portfolio;
                                    controller_delta = (0.0, 0.0, 0.0),
                                    load_delta = (0.0, 0.0),
                                    ibr_delta = 0.0,
                                    branch_delta = zeros(Float64, 46),
                                    bounds = :primary)
    control_limit, load_limit, ibr_limit, network_limit = if bounds == :primary
        (PRIMARY_CONTROL_BOUNDS, PRIMARY_LOAD_BOUNDS, PRIMARY_IBR_BOUNDS, PRIMARY_NETWORK_BOUNDS)
    elseif bounds == :discovery
        (0.20, PRIMARY_LOAD_BOUNDS, PRIMARY_IBR_BOUNDS, 0.20)
    elseif bounds == :secondary
        (0.25, 0.05, 0.05, 0.25)
    elseif bounds == :repair
        # A tuned baseline combined multiplicatively with the discovery
        # +/-20% screen can reach +/-32%; this is not a post-hoc expansion
        # of the repair coordinate itself.
        (0.35, 0.05, 0.05, 0.30)
    else
        throw(ArgumentError("bounds must be :primary, :discovery, :secondary, or :repair"))
    end
    c = Float64.(collect(controller_delta))
    length(c) == 3 || throw(ArgumentError("three control coordinates required"))
    all(abs.(c) .<= control_limit + 1e-12) ||
        throw(ArgumentError("controller coordinate outside primary box"))
    lp, lq = Float64.(collect(load_delta))
    abs(lp) <= load_limit + 1e-12 || throw(ArgumentError("load P coordinate outside selected box"))
    abs(lq) <= load_limit + 1e-12 || throw(ArgumentError("load Q coordinate outside selected box"))
    abs(ibr_delta) <= ibr_limit + 1e-12 || throw(ArgumentError("IBR coordinate outside selected box"))
    length(branch_delta) == 46 || throw(ArgumentError("46 branch coordinates required"))
    all(abs.(branch_delta) .<= network_limit + 1e-12) || throw(ArgumentError("branch coordinate outside selected box"))

    base = baseline_network()
    template = simple_gfldc_template(
        pll_scale = 1 + c[1],
        filter_scale = 1 + c[2],
        current_control_scale = 1 + c[3],
    )
    nw = replace_buses(base, portfolio; template = template)
    vertices, edges = copy_network_components(nw)
    _set_load_operating_point!(vertices, lp, lq)
    _set_replaced_dispatch!(vertices, portfolio, ibr_delta)
    apply_line_gamma(edges, branch_delta)
    out = Network(vertices, edges)
    set_jac_prototype!(out)
    return out
end

function _eq_residual(state)
    nw = extract_nw(state)
    du = zeros(Float64, length(uflat(state)))
    nw(du, uflat(state), pflat(state), state.t)
    return maximum(abs, du)
end

function _conditioning(sys, red)
    a = red.A
    sv = svdvals(a)
    a_cond = isempty(sv) ? NaN : maximum(sv) / max(minimum(sv), eps(Float64))
    gz_cond = NaN
    gz_min_sv = NaN
    if !(sys.M isa UniformScaling)
        md = diag(sys.M)
        cidx = findall(==(0), md)
        if !isempty(cidx)
            gz = sys.A[cidx, cidx]
            gsv = svdvals(gz)
            gz_min_sv = minimum(gsv)
            gz_cond = maximum(gsv) / max(gz_min_sv, eps(Float64))
        end
    end
    return (jacobian_condition = a_cond, g_z_condition = gz_cond,
            smallest_singular_value = isempty(sv) ? NaN : minimum(sv),
            g_z_smallest_singular_value = gz_min_sv)
end

function _critical_mode(state)
    sys = linearize_network(state)
    red = reduce_dae(sys)
    ee = eigen(red.A)
    λ = ComplexF64.(ee.values)
    gauge = findall(abs.(λ) .<= 1e-8)
    non = findall(abs.(λ) .> 1e-8)
    isempty(non) && throw(ArgumentError("no non-gauge mode"))
    k = non[argmax(real.(λ[non]))]
    V = ee.vectors
    W = inv(V)
    v = V[:, k]
    w = W[k, :]
    pf = abs.(w .* v)
    s = sum(pf)
    s > 0 && (pf ./= s)
    order = sortperm(pf; rev = true)
    top = order[1:min(8, length(order))]
    labels = isnothing(red.sym) ? String[] : string.(red.sym[top])
    converter_share = sum(pf[t] for t in top if occursin(r"SimpleGFLDC|PLL|CC1|gfl|dc", string(red.sym[t])); init = 0.0)
    machine_share = sum(pf[t] for t in top if occursin(r"machine|gov|avr|ω|delta|δ", string(red.sym[t])); init = 0.0)
    fam = if converter_share >= 0.25 && converter_share >= machine_share
        "converter-control oscillatory"
    elseif machine_share >= 0.25
        "electromechanical/control"
    elseif abs(imag(λ[k])) <= 1e-8
        "real aperiodic"
    else
        "slow/control unresolved"
    end
    condv = opnorm(V) * opnorm(W)
    return (sys = sys, red = red, eigenvalues = λ, critical_index = k,
            critical_eigenvalue = λ[k], critical_frequency_hz = abs(imag(λ[k])) / (2pi),
            damping_ratio = abs(real(λ[k])) / max(abs(λ[k]), eps(Float64)),
            gauge_count = length(gauge), nontrivial_count = length(non),
            critical_mode_family = fam, critical_state_labels = join(labels, "|"),
            critical_participation = join([string(red.sym[t], ":", pf[t]) for t in top], "|"),
            eigenvector_condition = condv, eigenvector = v)
end

"Return the complete confirmatory scalar/modal audit for a qualified state." 
function modal_report(state)
    cm = _critical_mode(state)
    cond = _conditioning(cm.sys, cm.red)
    return merge(cm, cond, (equilibrium_residual = _eq_residual(state),))
end

"Project the critical mode onto the common 39-bus complex voltage observable." 
function common_mode_shape(state)
    cm = _critical_mode(state)
    obs = Any[]
    for bus in 1:39
        push!(obs, VIndex(bus, :busbar₊u_r))
        push!(obs, VIndex(bus, :busbar₊u_i))
    end
    ps = NetworkDynamics.SII.parameter_symbols(state)
    isempty(ps) && throw(ArgumentError("no input channel available for common observable projection"))
    full = linearize_network(state; in = ps[1], out = obs)
    red = reduce_dae(full)
    y = ComplexF64.(red.C * cm.eigenvector)
    n = norm(y)
    n > 0 && (y ./= n)
    return y
end

function modal_assurance(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return NaN
    da = norm(a)
    db = norm(b)
    (da > 0 && db > 0) || return NaN
    return abs2(dot(conj(a), b)) / (da^2 * db^2)
end

struct ConfirmatoryCase
    status::String
    error_type::String
    error_message::String
    equilibrium_status::String
    state
    modal
end

function run_confirmatory_case(nw; label = "")
    connected = try
        is_connected(SimpleGraph(nw.im.g))
    catch
        false
    end
    connected || return ConfirmatoryCase("failed", "graph_disconnected", "network graph is disconnected",
                                         "graph_disconnected", nothing, nothing)
    local eq
    try
        eq = initialize_equilibrium(nw; sparse = false)
    catch err
        return ConfirmatoryCase("failed", "powerflow_or_initialization_exception", sprint(showerror, err),
                               "equilibrium_failed", nothing, nothing)
    end
    eq.powerflow_finite || return ConfirmatoryCase("failed", "powerflow_nonfinite", "PF state is nonfinite",
                                                   "equilibrium_failed", nothing, nothing)
    eq.state_finite || return ConfirmatoryCase("failed", "dynamic_initialization_nonfinite", "dynamic state is nonfinite",
                                                "equilibrium_failed", nothing, nothing)
    eq.fixed_point || return ConfirmatoryCase("failed", "fixed_point_failure", "fixed-point residual exceeds tolerance",
                                               "equilibrium_failed", eq.state, nothing)
    local modal
    try
        modal = modal_report(eq.state)
    catch err
        return ConfirmatoryCase("failed", "spectrum_failure", sprint(showerror, err),
                               "spectrum_failed", eq.state, nothing)
    end
    return ConfirmatoryCase("ok", "", "", "ok", eq.state, modal)
end

"Fast scalar holdout evaluator: same equilibrium and stability gates, no modal decomposition." 
function run_margin_case(nw)
    connected = try
        is_connected(SimpleGraph(nw.im.g))
    catch
        false
    end
    connected || return (status = "failed", error_type = "graph_disconnected",
        error_message = "network graph is disconnected", equilibrium_status = "graph_disconnected",
        alpha = NaN, margin = NaN, stable = false, equilibrium_residual = NaN)
    local eq
    try
        eq = initialize_equilibrium(nw; sparse = false)
        eq.powerflow_finite || return (status = "failed", error_type = "powerflow_nonfinite",
            error_message = "PF state is nonfinite", equilibrium_status = "equilibrium_failed",
            alpha = NaN, margin = NaN, stable = false, equilibrium_residual = NaN)
        eq.state_finite || return (status = "failed", error_type = "dynamic_initialization_nonfinite",
            error_message = "dynamic state is nonfinite", equilibrium_status = "equilibrium_failed",
            alpha = NaN, margin = NaN, stable = false, equilibrium_residual = NaN)
        eq.fixed_point || return (status = "failed", error_type = "fixed_point_failure",
            error_message = "fixed-point residual exceeds tolerance", equilibrium_status = "equilibrium_failed",
            alpha = NaN, margin = NaN, stable = false, equilibrium_residual = _eq_residual(eq.state))
        a = stability_audit(eq.state)
        return (status = "ok", error_type = "", error_message = "", equilibrium_status = "ok",
            alpha = a.max_real, margin = a.dynamic_margin, stable = a.stable,
            equilibrium_residual = _eq_residual(eq.state))
    catch err
        return (status = "failed", error_type = "powerflow_or_spectrum_exception",
            error_message = sprint(showerror, err), equilibrium_status = "equilibrium_failed",
            alpha = NaN, margin = NaN, stable = false, equilibrium_residual = NaN)
    end
end

function largest_load_bus()
    data = ieee39_data()
    ranks = sort([(abs(Float64(r.Pset)), Int(r.bus)) for r in eachrow(data.load)], by = x -> (-x[1], x[2]))
    return last(first(ranks))
end

"Attach the preregistered exact temporary constant-power-factor load pulse." 
function build_tds_network(portfolio; pulse = 0.01, load_bus = TDS_LOAD_BUS,
                            controller_delta = (0.0, 0.0, 0.0),
                            branch_delta = zeros(Float64, 46))
    nw = build_confirmatory_network(portfolio; controller_delta = controller_delta,
        branch_delta = branch_delta)
    vertices, edges = copy_network_components(nw)
    defaults = get_defaults_dict(vertices[load_bus])
    ps = only([s for s in keys(defaults) if occursin("Pset", string(s))])
    qs = only([s for s in keys(defaults) if occursin("Qset", string(s))])
    p0, q0 = defaults[ps], defaults[qs]
    affect = (u, p, ctx) -> begin
        factor = ctx.t < 1.05 ? 1 + pulse : 1.0
        p[ps] = p0 * factor
        p[qs] = q0 * factor
    end
    cb = PresetTimeComponentCallback([1.0, 1.1], ComponentAffect(affect, (), (ps, qs)))
    set_callback!(vertices[load_bus], cb)
    out = Network(vertices, edges)
    set_jac_prototype!(out)
    return out
end

function _try_state_value(s, bus, syms)
    for sym in syms
        try
            return Float64(s[VIndex(bus, sym)])
        catch
        end
    end
    return NaN
end

function _voltage_value(s, bus)
    ur = _try_state_value(s, bus, (:busbar₊u_r,))
    ui = _try_state_value(s, bus, (:busbar₊u_i,))
    return (isfinite(ur) && isfinite(ui)) ? hypot(ur, ui) : NaN
end

function _frequency_values(s)
    values = Float64[]
    for bus in 1:39
        ω = _try_state_value(s, bus, (:ctrld_gen₊machine₊ω, :machine₊ω,
                                      :gfl₊pll₊ω))
        isfinite(ω) && push!(values, ω * BASE_FREQ)
    end
    return values
end

function _settling_time(t, y, start; fraction = 0.02)
    idx = findall(>=(start), t)
    isempty(idx) && return NaN
    peak = maximum(abs.(y[idx]))
    peak <= eps(Float64) && return 0.0
    limit = fraction * peak
    for i in idx
        all(abs.(y[i:end]) .<= limit) && return t[i]
    end
    return NaN
end

function _log_rate(t, y; start = 1.1)
    idx = findall(i -> t[i] >= start && isfinite(y[i]) && y[i] > 1e-10, eachindex(t))
    length(idx) < 3 && return NaN
    x = Float64.(t[idx]); z = log.(Float64.(y[idx]))
    xm, zm = mean(x), mean(z)
    return sum((x .- xm) .* (z .- zm)) / max(sum((x .- xm).^2), eps(Float64))
end

"Run the preregistered nonlinear TDS and return the frozen observables." 
function simulate_tds_case(portfolio; pulse = 0.01, controller_delta = (0.0, 0.0, 0.0),
                           branch_delta = zeros(Float64, 46),
                           tspan = (0.0, 20.0), saveat = 0.01)
    nw = build_tds_network(portfolio; pulse = pulse, controller_delta = controller_delta,
        branch_delta = branch_delta)
    eq = initialize_equilibrium(nw; sparse = false)
    eq.powerflow_finite && eq.state_finite && eq.fixed_point ||
        throw(ArgumentError("TDS initial equilibrium is not qualified"))
    sref = eq.state
    vref = _voltage_value(sref, TDS_LOAD_BUS)
    prob = SciMLBase.ODEProblem(nw, sref, tspan)
    sol = SciMLBase.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5P();
        callback = get_callbacks(nw), initializealg = SciMLBase.NoInit(),
        saveat = saveat, abstol = 1e-8, reltol = 1e-8)
    t = Float64.(sol.t)
    states = [NetworkDynamics.NWState(sol, ti) for ti in t]
    v = [_voltage_value(s, TDS_LOAD_BUS) for s in states]
    spread = Float64[]
    for s in states
        f = _frequency_values(s)
        push!(spread, isempty(f) ? NaN : maximum(f) - minimum(f))
    end
    vdev = abs.(v .- vref)
    freq_finite = filter(isfinite, spread)
    v_finite = filter(isfinite, v)
    combined = [max(isfinite(spread[i]) ? spread[i] : 0.0,
                    isfinite(vdev[i]) ? BASE_FREQ * vdev[i] : 0.0) for i in eachindex(t)]
    return (solution = sol, t = t, voltage = v, frequency_spread_hz = spread,
        voltage_deviation_pu = vdev, voltage_min_pu = isempty(v_finite) ? NaN : minimum(v_finite),
        voltage_max_pu = isempty(v_finite) ? NaN : maximum(v_finite),
        max_frequency_spread_hz = isempty(freq_finite) ? NaN : maximum(freq_finite),
        max_voltage_deviation_pu = isempty(v_finite) ? NaN : maximum(vdev[isfinite.(v)]),
        frequency_settling_s = _settling_time(t, spread, 1.1),
        voltage_settling_s = _settling_time(t, vdev, 1.1),
        estimated_rate_s_inv = _log_rate(t, combined),
        tds_success = true)
end

"Compact row for the blocker/modal and holdout tables." 
function diagnostic_report(case::ConfirmatoryCase; case_id = "")
    if case.modal === nothing
        return (case_id = case_id, status = case.status, error_type = case.error_type,
                error_message = case.error_message, equilibrium_status = case.equilibrium_status,
                alpha = NaN, margin = NaN, stable = false, equilibrium_residual = NaN,
                critical_eigenvalue = "", critical_frequency_hz = NaN, damping_ratio = NaN,
                critical_mode_family = "", critical_state_labels = "", critical_participation = "",
                jacobian_condition = NaN, g_z_condition = NaN,
                smallest_singular_value = NaN, g_z_smallest_singular_value = NaN,
                eigenvector_condition = NaN)
    end
    m = case.modal
    alpha = real(m.critical_eigenvalue)
    return (case_id = case_id, status = case.status, error_type = case.error_type,
            error_message = case.error_message, equilibrium_status = case.equilibrium_status,
            alpha = alpha, margin = -alpha, stable = alpha < 0,
            equilibrium_residual = m.equilibrium_residual,
            critical_eigenvalue = string(m.critical_eigenvalue),
            critical_frequency_hz = m.critical_frequency_hz, damping_ratio = m.damping_ratio,
            critical_mode_family = m.critical_mode_family,
            critical_state_labels = m.critical_state_labels,
            critical_participation = m.critical_participation,
            jacobian_condition = m.jacobian_condition, g_z_condition = m.g_z_condition,
            smallest_singular_value = m.smallest_singular_value,
            g_z_smallest_singular_value = m.g_z_smallest_singular_value,
            eigenvector_condition = m.eigenvector_condition)
end
