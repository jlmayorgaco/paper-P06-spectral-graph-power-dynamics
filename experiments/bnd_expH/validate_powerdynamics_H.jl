using CSV, DataFrames, LinearAlgebra, TOML, SHA, Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPORT = joinpath(ROOT, "reports", "experiment_H")
const TABLES = joinpath(REPORT, "tables")
const CANDIDATE_PATH = joinpath(REPORT, "Z_H_FINAL.toml")
const SIGMA_REQUIRED = 0.05
const FREQ_HZ = 60.0
const EVENT_BUS = 16

# Fail closed before loading PowerDynamics: post-freeze validation is forbidden
# from changing either the candidate or its design-time provenance.
isfile(CANDIDATE_PATH) && isfile(CANDIDATE_PATH * ".sha256") || error("frozen ExpH candidate or SHA sidecar is missing")
frozen_sha = bytes2hex(sha256(read(CANDIDATE_PATH)))
strip(read(CANDIDATE_PATH * ".sha256", String)) == frozen_sha || error("frozen ExpH candidate SHA-256 mismatch")
candidate = TOML.parsefile(CANDIDATE_PATH)
candidate["candidate_frozen"] == true || error("candidate is not marked frozen")
candidate["design_used_PowerDynamics"] == false || error("PowerDynamics was used during design")
candidate["status"] == "FROZEN_ANALYTICAL_ROBUST_CANDIDATE" || error("unexpected candidate status")

# PowerDynamics is loaded only after the immutable candidate hash gate passes.
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase
include(joinpath(ROOT, "src", "pd39", "PD39.jl"))
using .PD39

function candidate_network(c, kp0, ki0)
    nw = PD39.baseline_network()
    rows = sort(c["generator"], by = x -> Int(x["bus"]))
    gain_mismatch = 0.0
    for row in rows
        bus = Int(row["bus"])
        rho = Float64(row["rho"])
        kp_scale = Float64(row["Kp"]) / kp0
        ki_scale = Float64(row["Ki"]) / ki0
        # PD39's frozen cross-fade uses one PLL frequency scale, giving Kp~b
        # and Ki~b^2. Record and gate the mismatch rather than silently
        # substituting an unrelated controller.
        gain_mismatch = max(gain_mismatch, abs(ki_scale - kp_scale^2))
        if rho == 1.0
            nw = PD39.replace_bus(nw, bus; template = PD39.simple_gfldc_template(pll_scale = kp_scale))
        elseif 0.0 < rho < 1.0
            nw = PD39.PD39Model.weighted_replacement_network(nw, bus, rho, kp_scale)
        elseif rho != 0.0
            error("invalid conversion fraction at bus $bus")
        end
    end
    gain_mismatch <= 1e-6 || error("candidate Kp/Ki cannot be represented by the validated PD39 one-scale PLL model; mismatch=$gain_mismatch")
    return nw, gain_mismatch
end

function frequency_state(state, bus)
    for sym in (:ctrld_gen₊machine₊ω, :machine₊ω, :gfl₊pll₊ω)
        try
            return Float64(state[VIndex(bus, sym)])
        catch
        end
    end
    return NaN
end

function voltage_pu(state, bus)
    try
        ur = Float64(state[VIndex(bus, :busbar₊u_r)])
        ui = Float64(state[VIndex(bus, :busbar₊u_i)])
        return hypot(ur, ui)
    catch
        return NaN
    end
end

function make_load_pulse_network(nw, bus, pulse)
    vertices, edges = PD39.PD39Model.copy_network_components(nw)
    defaults = get_defaults_dict(vertices[bus])
    pkeys = [s for s in keys(defaults) if occursin("Pset", string(s))]
    qkeys = [s for s in keys(defaults) if occursin("Qset", string(s))]
    length(pkeys) == 1 && length(qkeys) == 1 || error("event bus $bus does not expose one Pset/Qset parameter pair")
    ps, qs = only(pkeys), only(qkeys)
    p0, q0 = defaults[ps], defaults[qs]
    affect = (u, p, ctx) -> begin
        factor = ctx.t < 1.1 ? 1 + pulse : 1.0
        p[ps] = p0 * factor
        p[qs] = q0 * factor
    end
    callback = PresetTimeComponentCallback([1.0, 1.1], ComponentAffect(affect, (), (ps, qs)))
    set_callback!(vertices[bus], callback)
    out = Network(vertices, edges)
    set_jac_prototype!(out)
    return out
end

function run_tds(nw, eq, c, pulse; tspan = (0.0, 12.0), saveat = 0.02)
    pulse_nw = make_load_pulse_network(nw, EVENT_BUS, pulse)
    prob = SciMLBase.ODEProblem(pulse_nw, eq.state, tspan)
    sol = SciMLBase.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5P();
        callback = get_callbacks(pulse_nw), initializealg = SciMLBase.NoInit(),
        saveat = saveat, abstol = 1e-8, reltol = 1e-8)
    times = Float64.(sol.t)
    keep = [g for g in c["generator"] if Float64(g["epsilon"]) > 1e-8]
    weights = [Float64(g["H_MVA_s"]) * Float64(g["epsilon"]) for g in keep]
    refs = [frequency_state(eq.state, Int(g["bus"])) for g in keep]
    all(isfinite, refs) && sum(weights) > 0 || error("retained SG frequency states are unavailable")
    coi = Float64[]
    vdev = Float64[]
    vref = [voltage_pu(eq.state, b) for b in 1:39]
    for time in times
        state = NetworkDynamics.NWState(sol, time)
        omega = [frequency_state(state, Int(g["bus"])) for g in keep]
        all(isfinite, omega) || error("nonfinite retained SG frequency in TDS")
        push!(coi, FREQ_HZ * dot(weights, omega .- refs) / sum(weights))
        v = [voltage_pu(state, b) for b in 1:39]
        all(isfinite, v) || error("nonfinite bus voltage in TDS")
        push!(vdev, maximum(abs.(v .- vref)))
    end
    dt = diff(times)
    dcoi = diff(coi)
    ids = findall((dt .> 0) .& isfinite.(dt) .& isfinite.(dcoi))
    rocof = dcoi[ids] ./ dt[ids]
    isempty(rocof) && error("no finite positive-duration TDS intervals")
    all_voltage = [voltage_pu(NetworkDynamics.NWState(sol, t), b) for t in times for b in 1:39]
    return (; times, coi, voltage_deviation = vdev,
        max_abs_frequency_deviation_Hz = maximum(abs.(coi)),
        max_abs_rocof_Hz_s = maximum(abs.(rocof)),
        max_voltage_deviation_pu = maximum(vdev),
        voltage_min_pu = minimum(all_voltage), voltage_max_pu = maximum(all_voltage),
        successful = SciMLBase.successful_retcode(sol.retcode), rocof)
end

function main()
    mkpath(TABLES)
    domain = TOML.parsefile(joinpath(ROOT, "experiments", "bnd_expE", "configs", "DESIGN_DOMAIN_FROZEN.toml"))
    kp0 = Float64(domain["Kp_nom"])
    ki0 = Float64(domain["Ki_nom"])
    println("EXP_H_PD: candidate SHA-256 verified: ", frozen_sha); flush(stdout)
    baseline_eq = nothing
    candidate_eq = nothing
    baseline_audit = nothing
    candidate_audit = nothing
    candidate_nw = nothing
    gain_mismatch = NaN
    pd_error = ""
    try
        println("EXP_H_PD: initialize all-SG reference"); flush(stdout)
        baseline_nw = PD39.baseline_network()
        baseline_eq = PD39.initialize_equilibrium(baseline_nw; sparse = false, check = :error)
        baseline_audit = PD39.stability_audit(baseline_eq.state)
        println("EXP_H_PD: initialize frozen mixed candidate"); flush(stdout)
        candidate_nw, gain_mismatch = candidate_network(candidate, kp0, ki0)
        candidate_eq = PD39.initialize_equilibrium(candidate_nw; sparse = false, check = :error)
        candidate_audit = PD39.stability_audit(candidate_eq.state)
    catch err
        pd_error = sprint(showerror, err)
    end

    analytical_alpha = Float64(candidate["spectral_abscissa_s_inv"])
    candidate_alpha = candidate_audit === nothing ? NaN : candidate_audit.max_real
    eq_residual = candidate_eq === nothing ? NaN : PD39._eq_residual(candidate_eq.state)
    pd_status = candidate_eq !== nothing && candidate_eq.fixed_point && candidate_audit.finite ? "EVALUATED" : "BLOCKED_EXCEPTION"
    pd_ok = pd_status == "EVALUATED" && candidate_alpha <= -SIGMA_REQUIRED && abs(candidate_alpha - analytical_alpha) <= 1e-3
    pd_rows = DataFrame(
        candidate_sha256 = [frozen_sha], status = [pd_status],
        baseline_fixed_point = [baseline_eq === nothing ? false : baseline_eq.fixed_point],
        candidate_fixed_point = [candidate_eq === nothing ? false : candidate_eq.fixed_point],
        candidate_equilibrium_residual_inf = [eq_residual],
        candidate_raw_finite_poles = [candidate_audit === nothing ? 0 : length(candidate_audit.eigenvalues)],
        candidate_numerical_gauge_poles = [candidate_audit === nothing ? 0 : length(candidate_audit.gauge_eigenvalues)],
        alpha_analytical_s_inv = [analytical_alpha], alpha_PD_s_inv = [candidate_alpha],
        delta_alpha_PD_minus_analytical_s_inv = [candidate_alpha - analytical_alpha],
        strict_margin_required_s_inv = [-SIGMA_REQUIRED], meets_strict_margin = [candidate_alpha <= -SIGMA_REQUIRED],
        gain_scale_compatibility_error = [gain_mismatch], error = [pd_error])
    CSV.write(joinpath(TABLES, "TABLE_H20_powerdynamics_validation.csv"), pd_rows)
    if candidate_audit !== nothing
        CSV.write(joinpath(TABLES, "TABLE_H20_powerdynamics_candidate_poles.csv"), DataFrame(
            index = eachindex(candidate_audit.eigenvalues), real_s_inv = real.(candidate_audit.eigenvalues),
            imag_s_inv = imag.(candidate_audit.eigenvalues), abs_s_inv = abs.(candidate_audit.eigenvalues),
            is_numerical_gauge = abs.(candidate_audit.eigenvalues) .<= 1e-8))
    end

    tds_rows = NamedTuple[]
    traces = DataFrame[]
    tds_error = ""
    if candidate_eq !== nothing && candidate_eq.fixed_point && candidate_nw !== nothing
        for pulse in (5e-4, 1e-3)
            try
                println("EXP_H_PD: nonlinear load pulse ", pulse, " at bus ", EVENT_BUS); flush(stdout)
                result = run_tds(candidate_nw, candidate_eq, candidate, pulse)
                push!(tds_rows, (pulse_fraction = pulse, event_bus = EVENT_BUS,
                    status = result.successful ? "EVALUATED" : "SOLVER_FAILURE",
                    max_abs_frequency_deviation_Hz = result.max_abs_frequency_deviation_Hz,
                    max_abs_rocof_Hz_s = result.max_abs_rocof_Hz_s,
                    max_voltage_deviation_pu = result.max_voltage_deviation_pu,
                    voltage_min_pu = result.voltage_min_pu, voltage_max_pu = result.voltage_max_pu,
                    finite = all(isfinite, result.coi) && all(isfinite, result.rocof) && all(isfinite, result.voltage_deviation),
                    candidate_sha256 = frozen_sha, error = ""))
                push!(traces, DataFrame(time_s = result.times, pulse_fraction = fill(pulse, length(result.times)),
                    retained_sg_coi_frequency_deviation_Hz = result.coi,
                    max_voltage_deviation_pu = result.voltage_deviation))
            catch err
                msg = sprint(showerror, err)
                tds_error = isempty(tds_error) ? msg : tds_error * " | " * msg
                push!(tds_rows, (pulse_fraction = pulse, event_bus = EVENT_BUS,
                    status = "BLOCKED_EXCEPTION", max_abs_frequency_deviation_Hz = NaN,
                    max_abs_rocof_Hz_s = NaN, max_voltage_deviation_pu = NaN,
                    voltage_min_pu = NaN, voltage_max_pu = NaN, finite = false,
                    candidate_sha256 = frozen_sha, error = msg))
            end
        end
    else
        for pulse in (5e-4, 1e-3)
            push!(tds_rows, (pulse_fraction = pulse, event_bus = EVENT_BUS, status = "BLOCKED_NO_QUALIFIED_EQUILIBRIUM",
                max_abs_frequency_deviation_Hz = NaN, max_abs_rocof_Hz_s = NaN, max_voltage_deviation_pu = NaN,
                voltage_min_pu = NaN, voltage_max_pu = NaN, finite = false, candidate_sha256 = frozen_sha, error = pd_error))
        end
    end
    tds_df = DataFrame(tds_rows)
    if nrow(tds_df) == 2 && all(tds_df.finite)
        fsmall, flarge = tds_df.max_abs_frequency_deviation_Hz
        vsmall, vlarge = tds_df.max_voltage_deviation_pu
        tds_df[!, :relative_frequency_scaling_error] = fill(abs(flarge / max(fsmall, eps()) - 2.0) / 2.0, 2)
        tds_df[!, :relative_voltage_scaling_error] = fill(abs(vlarge / max(vsmall, eps()) - 2.0) / 2.0, 2)
    else
        tds_df[!, :relative_frequency_scaling_error] = fill(NaN, nrow(tds_df))
        tds_df[!, :relative_voltage_scaling_error] = fill(NaN, nrow(tds_df))
    end
    tds_ok = nrow(tds_df) == 2 && all(tds_df.status .== "EVALUATED") && all(tds_df.finite) &&
        maximum(tds_df.relative_frequency_scaling_error) <= 0.10 && maximum(tds_df.relative_voltage_scaling_error) <= 0.10
    CSV.write(joinpath(TABLES, "TABLE_H21_tds_validation.csv"), tds_df)
    isempty(traces) || CSV.write(joinpath(TABLES, "TABLE_H21_tds_trajectories.csv"), vcat(traces...))
    overall = pd_ok && tds_ok ? "PASS" : (pd_status == "EVALUATED" ? "PARTIAL" : "BLOCKED")
    open(joinpath(REPORT, "POSTFREEZE_VALIDATION_EXP_H.md"), "w") do io
        println(io, "# ExpH post-freeze PowerDynamics validation\n")
        println(io, "Candidate SHA-256: `$(frozen_sha)`. The SHA gate passed before PowerDynamics was loaded. The frozen candidate was not changed.\n")
        println(io, "Overall post-freeze status: **$(overall)**. Candidate equilibrium residual: $(eq_residual); PD spectral abscissa: $(candidate_alpha) s⁻¹; analytical abscissa: $(analytical_alpha) s⁻¹; difference: $(candidate_alpha - analytical_alpha) s⁻¹.\n")
        println(io, "Candidate has $(candidate_audit === nothing ? 0 : length(candidate_audit.eigenvalues)) finite raw poles and $(candidate_audit === nothing ? 0 : length(candidate_audit.gauge_eigenvalues)) numerical gauge poles. The comparison is between different model realizations; it is a margin-level check, not a pole-by-pole identity.\n")
        println(io, "Nonlinear load pulses at bus $(EVENT_BUS): 0.0005 and 0.001. TDS scaling gate: $(tds_ok ? "PASS" : "PARTIAL/BLOCKED"). See TABLE_H21_tds_validation.csv and TABLE_H21_tds_trajectories.csv.\n")
        println(io, "Controller compatibility error for the frozen Kp/Ki under PD39's shared PLL scale is $(gain_mismatch).\n")
        isempty(pd_error) || println(io, "PowerDynamics error: `$(pd_error)`\n")
        isempty(tds_error) || println(io, "TDS error: `$(tds_error)`\n")
        println(io, "This validation does not certify the analytical full-block uncertainty bound in the nonlinear PowerDynamics model, nor local/global optimality of the design.")
    end
    println("POWERDYNAMICS_VALIDATION: ", pd_status, " alpha=", candidate_alpha, " strict_margin=", candidate_alpha <= -SIGMA_REQUIRED)
    println("TDS_VALIDATION: ", tds_ok ? "PASS" : "PARTIAL/BLOCKED")
    println("EXP_H_POSTFREEZE_VALIDATION_STATUS: ", overall)
end

main()
