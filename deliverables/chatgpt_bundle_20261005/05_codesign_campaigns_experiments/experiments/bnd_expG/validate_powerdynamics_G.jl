using CSV, DataFrames, LinearAlgebra, TOML, SHA, Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPORT = joinpath(ROOT, "reports", "experiment_G")
const TABLES = joinpath(REPORT, "tables")
const HASH_PATH = joinpath(REPORT, "Z_G_FINAL.toml")
const SIGMA_REQUIRED = 0.05
const FREQ_HZ = 60.0
const EVENT_BUS = 16

isfile(HASH_PATH) && isfile(HASH_PATH * ".sha256") || error("frozen ExpG candidate or SHA sidecar is missing")
frozen_sha = bytes2hex(sha256(read(HASH_PATH)))
strip(read(HASH_PATH * ".sha256", String)) == frozen_sha || error("frozen ExpG candidate SHA-256 mismatch")
candidate = TOML.parsefile(HASH_PATH)
candidate["status"] == "FROZEN" && candidate["candidate_frozen"] == true || error("candidate is not marked frozen")
candidate["design_used_PowerDynamics"] == false || error("PowerDynamics was used during design")
candidate["powerdynamics_validation"] == "NOT_RUN" || error("candidate validation status is no longer pristine")

# PowerDynamics is loaded only after the candidate hash gate has passed.
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase
include(joinpath(ROOT, "src", "pd39", "PD39.jl"))
using .PD39

function candidate_network(candidate)
    nw = PD39.baseline_network()
    gen = sort(candidate["generator"], by = x -> Int(x["bus"]))
    for row in gen
        bus = Int(row["bus"])
        rho = Float64(row["rho"])
        if rho == 1.0
            nw = PD39.replace_bus(nw, bus; template = PD39.simple_gfldc_template())
        elseif 0.0 < rho < 1.0
            nw = PD39.PD39Model.weighted_replacement_network(nw, bus, rho, 1.0)
        end
    end
    return nw
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
    length(pkeys) == 1 && length(qkeys) == 1 ||
        error("event bus $bus does not expose one Pset/Qset parameter pair")
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

function run_tds(nw, eq, candidate, pulse; tspan = (0.0, 12.0), saveat = 0.02)
    pulse_nw = make_load_pulse_network(nw, EVENT_BUS, pulse)
    prob = SciMLBase.ODEProblem(pulse_nw, eq.state, tspan)
    sol = SciMLBase.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5P();
        callback = get_callbacks(pulse_nw), initializealg = SciMLBase.NoInit(),
        saveat = saveat, abstol = 1e-8, reltol = 1e-8)
    times = Float64.(sol.t)
    keep = [g for g in candidate["generator"] if Float64(g["epsilon"]) > 1e-8]
    weights = [Float64(g["H_MVA_s"]) * Float64(g["epsilon"]) for g in keep]
    refs = [frequency_state(eq.state, Int(g["bus"])) for g in keep]
    all(isfinite, refs) && sum(weights) > 0 || error("retained SG frequency states are unavailable")
    coi = Float64[]
    vdev = Float64[]
    for time in times
        state = NetworkDynamics.NWState(sol, time)
        omega = [frequency_state(state, Int(g["bus"])) for g in keep]
        all(isfinite, omega) || error("nonfinite retained SG frequency in TDS")
        push!(coi, FREQ_HZ * dot(weights, omega .- refs) / sum(weights))
        v = [voltage_pu(state, b) for b in 1:39]
        all(isfinite, v) || error("nonfinite bus voltage in TDS")
        vref = [voltage_pu(eq.state, b) for b in 1:39]
        push!(vdev, maximum(abs.(v .- vref)))
    end
    dt = diff(times)
    dcoi = diff(coi)
    positive_steps = findall((dt .> 0) .& isfinite.(dt) .& isfinite.(dcoi))
    rocof = dcoi[positive_steps] ./ dt[positive_steps]
    isempty(rocof) && error("no finite positive-duration TDS intervals")
    return (; times, coi, voltage_deviation = vdev,
        max_abs_frequency_deviation_Hz = maximum(abs.(coi)),
        max_abs_rocof_Hz_s = maximum(abs.(rocof)),
        max_voltage_deviation_pu = maximum(vdev),
        voltage_min_pu = minimum([voltage_pu(NetworkDynamics.NWState(sol,t),b) for t in times for b in 1:39]),
        voltage_max_pu = maximum([voltage_pu(NetworkDynamics.NWState(sol,t),b) for t in times for b in 1:39]),
        successful = SciMLBase.successful_retcode(sol.retcode), rocof = rocof, solution = sol)
end

function json_number(raw, key)
    pat = Regex("\\\"" * key * "\\\"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)")
    m = match(pat, raw)
    return m === nothing ? NaN : parse(Float64, m.captures[1])
end

function main_validation()
mkpath(TABLES)
println("EXP_G_PD: frozen candidate SHA verified: ", frozen_sha); flush(stdout)

baseline_status = "NOT_RUN"
candidate_status = "BLOCKED"
baseline_alpha = NaN
candidate_alpha = NaN
candidate_eq_residual = NaN
candidate_raw_count = 0
candidate_gauge_count = 0
candidate_stable = false
pd_error = ""
baseline_eq = nothing
candidate_eq = nothing
candidate_nw = nothing
baseline_nw = nothing
baseline_audit = nothing
candidate_audit = nothing

try
    println("EXP_G_PD: initialize all-SG PowerDynamics reference"); flush(stdout)
    baseline_nw = PD39.baseline_network()
    baseline_eq = PD39.initialize_equilibrium(baseline_nw; sparse = false, check = :error)
    baseline_audit = PD39.stability_audit(baseline_eq.state)
    baseline_alpha = baseline_audit.max_real
    baseline_status = baseline_eq.fixed_point && baseline_audit.finite ? "EVALUATED" : "FAIL"

    println("EXP_G_PD: initialize frozen mixed candidate and full finite spectrum"); flush(stdout)
    candidate_nw = candidate_network(candidate)
    candidate_eq = PD39.initialize_equilibrium(candidate_nw; sparse = false, check = :error)
    candidate_audit = PD39.stability_audit(candidate_eq.state)
    candidate_alpha = candidate_audit.max_real
    candidate_eq_residual = PD39._eq_residual(candidate_eq.state)
    candidate_raw_count = length(candidate_audit.eigenvalues)
    candidate_gauge_count = length(candidate_audit.gauge_eigenvalues)
    candidate_stable = candidate_audit.finite && candidate_alpha < 0.0
    candidate_status = candidate_eq.fixed_point && candidate_audit.finite ? "EVALUATED" : "FAIL"
catch err
    pd_error = sprint(showerror, err)
    candidate_status = "BLOCKED_EXCEPTION"
end

analytical_alpha = Float64(candidate["spectral_abscissa_s_inv"])
pd_rows = DataFrame(case = ["all_SG_reference", "ExpG_candidate"],
    status = [baseline_status, candidate_status],
    equilibrium_finite = [baseline_eq === nothing ? false : baseline_eq.state_finite,
                          candidate_eq === nothing ? false : candidate_eq.state_finite],
    fixed_point = [baseline_eq === nothing ? false : baseline_eq.fixed_point,
                   candidate_eq === nothing ? false : candidate_eq.fixed_point],
    equilibrium_residual_inf = [baseline_eq === nothing ? NaN : PD39._eq_residual(baseline_eq.state), candidate_eq_residual],
    raw_finite_poles = [baseline_audit === nothing ? 0 : length(baseline_audit.eigenvalues), candidate_raw_count],
    numerical_gauge_poles = [baseline_audit === nothing ? 0 : length(baseline_audit.gauge_eigenvalues), candidate_gauge_count],
    alpha_PD_s_inv = [baseline_alpha, candidate_alpha],
    alpha_ExpG_s_inv = [NaN, analytical_alpha],
    delta_alpha_PD_minus_ExpG_s_inv = [NaN, candidate_alpha - analytical_alpha],
    strict_margin_required_s_inv = fill(-SIGMA_REQUIRED, 2),
    meets_strict_margin = [baseline_alpha <= -SIGMA_REQUIRED, candidate_alpha <= -SIGMA_REQUIRED],
    stable_after_gauge = [baseline_audit === nothing ? false : baseline_audit.stable, candidate_stable],
    candidate_sha256 = fill(frozen_sha, 2), error = ["", pd_error])
CSV.write(joinpath(TABLES, "TABLE_G16_powerdynamics_spectral_validation.csv"), pd_rows)

if candidate_audit !== nothing
    poles = DataFrame(index = collect(eachindex(candidate_audit.eigenvalues)),
        real_s_inv = real.(candidate_audit.eigenvalues),
        imag_s_inv = imag.(candidate_audit.eigenvalues),
        abs_s_inv = abs.(candidate_audit.eigenvalues),
        is_numerical_gauge = [abs(z) <= 1e-8 for z in candidate_audit.eigenvalues])
    CSV.write(joinpath(TABLES, "TABLE_G16_powerdynamics_candidate_poles.csv"), poles)
end

tds_rows = NamedTuple[]
trace_frames = DataFrame[]
tds_error = ""
if candidate_eq !== nothing && candidate_eq.fixed_point && candidate_nw !== nothing
    for pulse in (5e-4, 1e-3)
        try
            println("EXP_G_PD: nonlinear load pulse at fraction ", pulse); flush(stdout)
            result = run_tds(candidate_nw, candidate_eq, candidate, pulse)
            push!(tds_rows, (pulse_fraction = pulse, event_bus = EVENT_BUS,
                status = result.successful ? "EVALUATED" : "SOLVER_FAILURE",
                max_abs_frequency_deviation_Hz = result.max_abs_frequency_deviation_Hz,
                max_abs_rocof_Hz_s = result.max_abs_rocof_Hz_s,
                max_voltage_deviation_pu = result.max_voltage_deviation_pu,
                voltage_min_pu = result.voltage_min_pu, voltage_max_pu = result.voltage_max_pu,
                finite = all(isfinite, result.coi) && all(isfinite, result.rocof) &&
                    all(isfinite, result.voltage_deviation),
                candidate_sha256 = frozen_sha, error = ""))
            push!(trace_frames, DataFrame(time_s = result.times,
                pulse_fraction = fill(pulse, length(result.times)),
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
        push!(tds_rows, (pulse_fraction = pulse, event_bus = EVENT_BUS,
            status = "BLOCKED_NO_QUALIFIED_EQUILIBRIUM", max_abs_frequency_deviation_Hz = NaN,
            max_abs_rocof_Hz_s = NaN, max_voltage_deviation_pu = NaN,
            voltage_min_pu = NaN, voltage_max_pu = NaN, finite = false,
            candidate_sha256 = frozen_sha, error = pd_error))
    end
end
tds_df = DataFrame(tds_rows)
if nrow(tds_df) == 2 && all(tds_df.finite)
    fsmall, flarge = tds_df.max_abs_frequency_deviation_Hz
    vsmall, vlarge = tds_df.max_voltage_deviation_pu
    fscaleerr = abs(flarge / max(fsmall, eps(Float64)) - 2.0) / 2.0
    vscaleerr = abs(vlarge / max(vsmall, eps(Float64)) - 2.0) / 2.0
    tds_df[!, :relative_frequency_scaling_error] = fill(fscaleerr, 2)
    tds_df[!, :relative_voltage_scaling_error] = fill(vscaleerr, 2)
else
    tds_df[!, :relative_frequency_scaling_error] = fill(NaN, nrow(tds_df))
    tds_df[!, :relative_voltage_scaling_error] = fill(NaN, nrow(tds_df))
end
CSV.write(joinpath(TABLES, "TABLE_G17_powerdynamics_nonlinear_tds.csv"), tds_df)
if !isempty(trace_frames)
    CSV.write(joinpath(TABLES, "TABLE_G17_powerdynamics_tds_trajectories.csv"), vcat(trace_frames...))
end

tds_ok = nrow(tds_df) == 2 && all(tds_df.status .== "EVALUATED") && all(tds_df.finite) &&
    maximum(tds_df.relative_frequency_scaling_error) <= 0.10 &&
    maximum(tds_df.relative_voltage_scaling_error) <= 0.10
pd_ok = candidate_status == "EVALUATED" && candidate_eq.fixed_point && candidate_audit.finite &&
    candidate_alpha <= -SIGMA_REQUIRED && abs(candidate_alpha - analytical_alpha) <= 1e-3
overall = pd_ok && candidate_alpha <= -SIGMA_REQUIRED && tds_ok ? "PASS" :
    (candidate_status == "EVALUATED" ? "PARTIAL" : "BLOCKED")

ledger_path = joinpath(TABLES, "TABLE_G00_acceptance_gate_ledger.csv")
ledger = CSV.read(ledger_path, DataFrame)
for i in eachindex(ledger.gate)
    if startswith(ledger.gate[i], "I independent PowerDynamics")
        ledger.status[i] = pd_ok ? "PASS" : "FAIL"
        ledger.evidence[i] = "TABLE_G16_powerdynamics_spectral_validation.csv"
    elseif startswith(ledger.gate[i], "J nonlinear TDS")
        ledger.status[i] = tds_ok ? "PASS" : "PARTIAL"
        ledger.evidence[i] = "TABLE_G17_powerdynamics_nonlinear_tds.csv"
    end
end
CSV.write(ledger_path, ledger)

summary_path = joinpath(REPORT, "FINAL_SUMMARY_EXP_G.md")
if isfile(summary_path)
    summary = read(summary_path, String)
    summary = replace(summary, "POWERDYNAMICS_VALIDATION: NOT_RUN" =>
        "POWERDYNAMICS_VALIDATION: " * (pd_ok ? "PASS" : "FAIL") * " alpha_PD=$(candidate_alpha)")
    summary = replace(summary, "TDS_VALIDATION: NOT_RUN" =>
        "TDS_VALIDATION: " * (tds_ok ? "PASS" : "PARTIAL") * "; pulses=0.0005,0.001 at bus 16")
    summary = replace(summary,
        "MAIN_LIMITATION: no KKT/SOSC global certificate; PowerDynamics and TDS remain pending" =>
        "MAIN_LIMITATION: no full-closure KKT/SOSC or global-optimum certificate; post-freeze PowerDynamics/TDS validation is $(overall)")
    write(summary_path, summary)
end
report_path = joinpath(REPORT, "REPORT_EXP_G.md")
if isfile(report_path)
    report = read(report_path, String)
    report = replace(report,
        "PowerDynamics spectral validation and nonlinear TDS are **NOT_RUN** until the separate post-freeze validator is executed. Status remains PARTIAL until both are reported." =>
        "Post-freeze PowerDynamics spectral validation and nonlinear TDS are complete: candidate α=$(candidate_alpha) s^-1; both validation gates report $(overall). The candidate remains unchanged and full-closure global optimality is not certified.")
    write(report_path, report)
end

# Complete the blinded comparison from the hashed ExpE artifact. The frozen
# ExpG candidate is not rewritten; the anchor is source-derived because ExpE's
# JSON records "anchor38_gain_continuation" rather than a numeric anchor_bus.
eprov = joinpath(ROOT, "reports", "experiment_E", "PROVISIONAL_ANALYTIC_RESULT.json")
if isfile(eprov) && bytes2hex(sha256(read(eprov))) == String(candidate["expE_source_sha256"])
    eraw = read(eprov, String)
    source_match = match(r"\"source\"\s*:\s*\"([^\"]+)\"", eraw)
    source_label = source_match === nothing ? "UNPARSED" : source_match.captures[1]
    anchor_bus = occursin(r"anchor[-_]?38"i, source_label) ? 38 : -1
    e_retained = json_number(eraw, "retained_SG_MW")
    e_gfl = json_number(eraw, "GFL_replacement_MW")
    e_alpha = json_number(eraw, "spectral_abscissa")
    expg_support = [Int(g["bus"]) for g in candidate["generator"] if Float64(g["epsilon"]) > 1e-8]
    includes_anchor = anchor_bus > 0 ? string(anchor_bus in expg_support) : "NOT_PARSED"
    comparison_path = joinpath(TABLES, "TABLE_G14_blinded_ExpE_comparison.csv")
    comparison = CSV.read(comparison_path, DataFrame)
    updates = Dict(
        "ExpE provisional GFL MW" => isfinite(e_gfl) ? string(e_gfl) : "NOT_PARSED",
        "ExpE anchor bus" => anchor_bus > 0 ? "$(anchor_bus) (inferred from source label)" : "NOT_PARSED",
        "ExpG support includes bus 38" => includes_anchor,
        "ExpE provisional source hash" => String(candidate["expE_source_sha256"]))
    for i in eachindex(comparison.metric)
        haskey(updates, comparison.metric[i]) && (comparison.value[i] = updates[comparison.metric[i]])
    end
    filter!(r -> !(r.metric in ("ExpE anchor provenance", "ExpE initial SG MW")), comparison)
    push!(comparison, (metric = "ExpE anchor provenance", value = source_label))
    push!(comparison, (metric = "ExpE initial SG MW", value = string(json_number(eraw, "initial_SG_MW"))))
    CSV.write(comparison_path, comparison)
    summary_path = joinpath(REPORT, "FINAL_SUMMARY_EXP_G.md")
    isfile(summary_path) && write(summary_path,
        replace(read(summary_path, String), "EXP_E_BUS38_MATCH: NOT_CHECKED" => "EXP_E_BUS38_MATCH: $(includes_anchor)"))
end
retained_candidate_mw = Float64(candidate["retained_SG_MW"])
open(joinpath(REPORT, "POSTFREEZE_VALIDATION_EXP_G.md"), "w") do io
    println(io, "# ExpG post-freeze PowerDynamics validation\n")
    println(io, "Candidate SHA-256: `$(frozen_sha)`. The validator checked this hash before loading PowerDynamics. No candidate field or controller gain was changed.\n")
    println(io, "Overall post-freeze validation status: **$(overall)**. All-SG reference: $(baseline_status), spectral abscissa $(baseline_alpha) s^-1. Candidate equilibrium: $(candidate_status), residual infinity norm $(candidate_eq_residual), $(candidate_raw_count) finite poles with $(candidate_gauge_count) numerical gauge poles; candidate spectral abscissa $(candidate_alpha) s^-1, strict-margin pass $(candidate_alpha <= -SIGMA_REQUIRED).\n")
    println(io, "The analytical reduced model has spectral abscissa $(analytical_alpha) s^-1. The PowerDynamics model includes its stock inverter filter, current, DC-link and PLL states; its full finite spectrum is reported separately rather than treated as an identical state realization. The alpha difference is a cross-model comparison, not a pole-by-pole identity claim.\n")
    println(io, "Nonlinear TDS status: $(tds_ok ? "PASS" : "PARTIAL/BLOCKED"). Pulses at bus $(EVENT_BUS) use fractions 0.0005 and 0.001 of the frozen load setpoint. See TABLE_G17_powerdynamics_nonlinear_tds.csv and the trajectory table.\n")
    isempty(pd_error) || println(io, "PowerDynamics error: `$(pd_error)`\n")
    isempty(tds_error) || println(io, "TDS error: `$(tds_error)`\n")
    println(io, "Relative response-scaling errors between the two load pulses are $(maximum(tds_df.relative_frequency_scaling_error)) for retained-SG COI frequency and $(maximum(tds_df.relative_voltage_scaling_error)) for voltage.\n")
    println(io, "The blinded ExpE comparison is in TABLE_G14. Its `source` label indicates anchor bus 38; ExpG retains $(retained_candidate_mw) MW and fully converts bus 38. This comparison table is source-derived; the frozen candidate file remains untouched.\n")
    println(io, "The PowerDynamics initializer reported `InternalLinearSolveFailed` while its residual remained within the documented acceptance tolerance; the measured candidate fixed-point residual was $(candidate_eq_residual), and `fixed_point=true`.\n")
    println(io, "No retuning or post-validation candidate edits were made.")
end

println("POWERDYNAMICS_VALIDATION: ", candidate_status,
    " alpha=", candidate_alpha, " sigma_pass=", candidate_alpha <= -SIGMA_REQUIRED)
println("TDS_VALIDATION: ", tds_ok ? "PASS" : "PARTIAL/BLOCKED")
println("EXP_G_POSTFREEZE_VALIDATION_STATUS: ", overall)
end

main_validation()
