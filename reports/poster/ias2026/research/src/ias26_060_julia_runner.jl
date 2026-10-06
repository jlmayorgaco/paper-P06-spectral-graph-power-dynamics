"""Run the frozen IAS26-060 scenario manifest through the exact custom Julia TX4 model.

The canonical TX4 source is loaded in memory without editing its frozen file.
Only the compile-time P4_G binding is made runtime-selectable so the two
preregistered treatments (g=0.03625 and g=0.25) can be evaluated. All device,
network, and DAE equations remain the unchanged canonical implementation.
"""

using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Printf
using SHA
using Statistics

function arg_value(flag::String; default=nothing)
    i = findfirst(==(flag), ARGS)
    i === nothing && return default
    i < length(ARGS) || error("$flag requires a value")
    ARGS[i + 1]
end

function sha256_file(path::String)
    bytes2hex(sha256(read(path)))
end

const REPO_ROOT = abspath(arg_value("--repo-root"; default=pwd()))
const TX4_SOURCE_PATH = joinpath(REPO_ROOT, "code", "tx4", "run_tx4_exact_p4_ieee39.jl")
const TX4_DATA_PATH = joinpath(REPO_ROOT, "code", "tx4", "frozen_ieee39_data.jl")
const TX4 = let
    source = read(TX4_SOURCE_PATH, String)
    needle = "const P4_G = 0.03625"
    count(line -> strip(line) == needle, split(source, '\n')) == 1 || error("expected one frozen P4_G binding in exact Julia runner")
    include_call = "include(joinpath(@__DIR__,"
    length(findall(include_call, source)) == 2 || error("unexpected include structure in exact Julia runner")
    adapted_source = replace(source, needle => "P4_G = 0.03625")
    adapted_source = replace(adapted_source, include_call => "Base.include(@__MODULE__, joinpath(@__DIR__,")
    model_module = Module(:TX4, true, true)
    Base.include_string(model_module, adapted_source, TX4_SOURCE_PATH)
    model_module
end

function parse_numeric_bus_map(raw::AbstractString)
    result = Dict{Int,Float64}()
    for m in eachmatch(r"\"([0-9]+)\"\s*:\s*(-?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)", raw)
        result[parse(Int, m.captures[1])] = parse(Float64, m.captures[2])
    end
    isempty(result) && error("could not parse numeric bus map")
    result
end

function parse_int_array(raw::AbstractString)
    [parse(Int, m.match) for m in eachmatch(r"[0-9]+", raw)]
end

function load_scenarios(path::String)
    [row for row in CSV.File(path; types=String, normalizenames=false)]
end

function select_smoke_rows(rows)
    scales = [parse(Float64, row[:load_scale_L_s]) for row in rows]
    indices = [argmin(scales), argmax(scales), argmin(abs.(scales .- 1.0)), 1, length(rows)]
    unique_indices = unique(indices)
    length(unique_indices) == 5 || error("declared five smoke selectors did not resolve to five unique IDs: $unique_indices")
    rows[unique_indices]
end

function load_scenario!(TX, row)
    p_load = parse_numeric_bus_map(row[:P_L_s_by_bus_pu_json])
    q_load = parse_numeric_bus_map(row[:Q_L_s_by_bus_pu_json])
    p_gen = parse_numeric_bus_map(row[:P_G_scheduled_by_bus_pu_json])
    load_buses = sort!(collect(keys(p_load)))
    load_buses == sort!(collect(keys(q_load))) || error("P/Q load bus ordering mismatch for $(row[:scenario_id])")

    empty!(TX.LOAD_ROWS)
    empty!(TX.LOAD_BY_BUS)
    for bus in load_buses
        p, q = p_load[bus], q_load[bus]
        push!(TX.LOAD_ROWS, (Float64(bus), p, q))
        TX.LOAD_BY_BUS[bus] = complex(p, q)
    end

    old_pv = Dict(Int(round(r[1])) => r for r in TX.PV_ROWS)
    length(p_gen) == length(old_pv) || error("generator bus-set mismatch for $(row[:scenario_id])")
    sort!(collect(keys(p_gen))) == sort!(collect(keys(old_pv))) || error("generator bus ordering mismatch for $(row[:scenario_id])")
    empty!(TX.PV_ROWS)
    for bus in sort!(collect(keys(old_pv)))
        r = old_pv[bus]
        push!(TX.PV_ROWS, (r[1], p_gen[bus], r[3], r[4], r[5], r[6], r[7], r[8]))
    end
    empty!(TX.PV_BY_BUS)
    for r in TX.PV_ROWS
        TX.PV_BY_BUS[Int(round(r[1]))] = r
    end
    TX.SCHEDULED .= TX.scheduled_power()
    (p_load=p_load, q_load=q_load, p_gen=p_gen)
end

function set_g!(TX, g::Float64)
    (g == 0.03625 || g == 0.25) || error("unregistered controller treatment g=$g")
    TX.P4_G = g
end

function physical_screen(TX, v, p_gen, settings)
    generation = v .* conj.(TX.YBUS * v)
    for (bus, load) in TX.LOAD_BY_BUS
        generation[bus] += load
    end
    tol = settings.limit_tolerance
    vmin, vmax = minimum(abs.(v)), maximum(abs.(v))
    failures = String[]
    vmin < settings.v_min - tol && push!(failures, "BUS_VOLTAGE_LIMIT")
    vmax > settings.v_max + tol && push!(failures, "BUS_VOLTAGE_LIMIT")
    for r in TX.PV_ROWS
        bus = Int(round(r[1]))
        p, q = real(generation[bus]), imag(generation[bus])
        sn_pu = TX.MACHINE_BY_BUS[bus][2] / TX.SYSTEM_BASE_MVA
        p < r[6] - tol || p > r[5] + tol ? push!(failures, "GENERATOR_P_LIMIT") : nothing
        q < r[8] - tol || q > r[7] + tol ? push!(failures, "GENERATOR_Q_LIMIT") : nothing
        hypot(p, q) > sn_pu + tol && push!(failures, "GENERATOR_MVA_LIMIT")
    end
    slack_bus = Int(round(TX.SLACK_ROW[1]))
    slack_p, slack_q = real(generation[slack_bus]), imag(generation[slack_bus])
    slack_p < settings.slack_p_min - tol || slack_p > settings.slack_p_max + tol ? push!(failures, "SLACK_P_LIMIT") : nothing
    slack_sn = Float64(TX.SLACK_ROW[4]) / TX.SYSTEM_BASE_MVA
    hypot(slack_p, slack_q) > slack_sn + tol && push!(failures, "SLACK_MVA_LIMIT")
    (failures=unique(failures), vmin=vmin, vmax=vmax,
     gen_pq=Dict(bus => (real(generation[bus]), imag(generation[bus])) for bus in sort!(collect(keys(p_gen)))),
     slack_p=slack_p, slack_q=slack_q)
end

function mode_families(values::AbstractVector{<:Complex})
    used = falses(length(values))
    families = NamedTuple[]
    for i in eachindex(values)
        used[i] && continue
        λ = values[i]
        partner = i
        if abs(imag(λ)) > 1e-10
            candidates = [j for j in eachindex(values) if j != i && !used[j]]
            isempty(candidates) || (partner = candidates[argmin(abs.(values[candidates] .- conj(λ)))])
        end
        used[i] = true
        used[partner] = true
        representative = imag(λ) >= 0 ? i : partner
        push!(families, (indices=sort(unique([i, partner])), representative=representative,
                         lambda=values[representative], frequency=abs(imag(values[representative]))/(2pi)))
    end
    sort!(families; by=f -> (f.frequency, real(f.lambda), imag(f.lambda)))
    families
end

function evaluate_case(TX, replaced::Set{Int}, g::Float64, pf)
    v, pf_residual, pf_iterations = pf
    set_g!(TX, g)
    d, x = TX.build_dae(v, replaced)
    z = zeros(Float64, 2 * length(TX.BUS_IDS))
    z[1:2:end] = real.(v)
    z[2:2:end] = imag.(v)
    f0, g0 = TX.f_residual(d, x, z), TX.g_residual(d, x, z)
    norm_f, norm_g = norm(f0, Inf), norm(g0, Inf)
    all(isfinite, x) && all(isfinite, z) && all(isfinite, f0) && all(isfinite, g0) || error("DAE_NONFINITE")

    t0 = time()
    fx, fz, gx, gz = TX.jacobians(d, x, z)
    gz_condition = cond(gz)
    isfinite(gz_condition) || error("EIGENSOLVER_FAILURE: nonfinite algebraic Jacobian condition")
    A = fx - fz * (gz \ gx)
    all(isfinite, A) || error("SPECTRUM_NONFINITE: reduced Jacobian contains nonfinite entries")
    eig = eigen(A)
    values, vectors = ComplexF64.(eig.values), ComplexF64.(eig.vectors)
    all(isfinite, values) && all(isfinite, vectors) || error("SPECTRUM_NONFINITE")
    transverse = findall(abs.(values) .>= TX.GAUGE_TOL)
    isempty(transverse) && error("EIGENSOLVER_FAILURE: empty transverse spectrum")
    critical_index = transverse[argmax(real.(values[transverse]))]
    critical = values[critical_index]
    families = mode_families(values)
    family_index = findfirst(f -> critical_index in f.indices, families)
    family = families[family_index]
    critical_mode_id = "FAMILY_$(lpad(string(family_index), 3, '0'))"

    target = Float64[]
    for i in transverse
        λ = values[i]
        hz = imag(λ)/(2pi)
        if imag(λ) >= -1e-10 && hz >= 0.3 - 1e-10 && hz <= 1.5 + 1e-10 && abs(λ) > 1e-3
            push!(target, real(λ))
        end
    end
    alpha_omega = isempty(target) ? NaN : maximum(target)
    a_norm_2 = opnorm(A, 2)
    eig_residuals = [norm(A * vectors[:, i] - values[i] * vectors[:, i], 2) /
                     (max(a_norm_2, 1e-300) * max(norm(vectors[:, i], 2), 1e-300)) for i in eachindex(values)]
    eig_elapsed = time() - t0
    alpha_perp = maximum(real.(values[transverse]))

    (status="NUMERICALLY_COMPLETE", failure_code="", g=g,
     state_count=length(x), transverse_dimension=length(transverse),
     pf_residual=pf_residual, pf_iterations=pf_iterations,
     dae_f_residual=norm_f, dae_g_residual=norm_g, gz_condition=gz_condition,
     alpha_perp=alpha_perp, alpha_omega=alpha_omega,
     critical_lambda=critical, critical_frequency_hz=abs(imag(critical))/(2pi),
     critical_damping_ratio=-real(critical)/max(abs(critical), 1e-300),
     critical_mode_id=critical_mode_id, critical_eigen_index=critical_index,
     eigenpair_residual_max=maximum(eig_residuals), eigen_elapsed_s=eig_elapsed,
     x=x, z=z, v=v, A=A, eigenvalues=values, eigenvectors=vectors,
     transverse_mask=BitVector([i in transverse for i in eachindex(values)]),
     eigenpair_residuals=eig_residuals, mode_family_indices=[f.indices for f in families],
     mode_family_eigenvalues=ComplexF64[f.lambda for f in families],
     mode_family_frequencies_hz=Float64[f.frequency for f in families])
end

function portfolio_spec()
    order = [30, 33, 35, 37]
    specs = [("BASE", Int[]), ("30", [30]), ("33", [33]), ("35", [35]), ("37", [37]),
             ("30+33", [30,33]), ("30+35", [30,35]), ("30+37", [30,37]),
             ("33+35", [33,35]), ("33+37", [33,37]), ("35+37", [35,37]),
             ("30+33+35", [30,33,35]), ("30+33+37", [30,33,37]),
             ("30+35+37", [30,35,37]), ("33+35+37", [33,35,37]), ("30+33+35+37", [30,33,35,37])]
    function mask_of(buses)
        mask = 0
        for bus in buses
            mask |= 1 << (findfirst(==(bus), order) - 1)
        end
        mask
    end
    [(id=id, buses=Set(buses), mask=mask_of(buses)) for (id,buses) in specs]
end

function base_case_hash(TX, result, scenario_id::String, portfolio_id::String, g::Float64)
    values = vcat(real.(result.v), imag.(result.v), result.x, result.z)
    payload = join((@sprintf("%.17g", x) for x in values), ",") * "|$scenario_id|$portfolio_id|$(@sprintf("%.17g", g))"
    bytes2hex(sha256(codeunits(payload)))
end

function write_scenario_spectrum(path::String, scenario_id::String, case_payloads)
    JLD2.jldopen(path, "w"; compress=true) do file
        file["metadata/scenario_id"] = scenario_id
        file["metadata/case_count"] = length(case_payloads)
        file["metadata/case_prefixes"] = ["cases/$(lpad(string(item.case_index), 2, '0'))_$(item.portfolio_id)_g$(replace(string(item.g), "."=>"p"))" for item in case_payloads]
        for (i, item) in enumerate(case_payloads)
            prefix = "cases/$(lpad(string(item.case_index), 2, '0'))_$(item.portfolio_id)_g$(replace(string(item.g), "."=>"p"))"
            file["$prefix/eigenvalues"] = item.eigenvalues
            file["$prefix/eigenvectors"] = item.eigenvectors
            file["$prefix/A_full"] = item.A
            file["$prefix/state_x"] = item.x
            file["$prefix/algebraic_z"] = item.z
            file["$prefix/bus_voltage"] = item.v
            file["$prefix/transverse_mask"] = item.transverse_mask
            file["$prefix/eigenpair_relative_residuals"] = item.eigenpair_residuals
            file["$prefix/mode_family_indices"] = item.mode_family_indices
            file["$prefix/mode_family_eigenvalues"] = item.mode_family_eigenvalues
        end
    end
end

function validate_scenario_spectrum(path::String, scenario_id::String, expected_cases::Int)
    JLD2.jldopen(path, "r") do file
        file["metadata/scenario_id"] == scenario_id || error("JLD2 scenario ID mismatch")
        file["metadata/case_count"] == expected_cases || error("JLD2 case count mismatch")
        if expected_cases > 0
            prefix = first(file["metadata/case_prefixes"])
            first_values = file["$prefix/eigenvalues"]
            first_vectors = file["$prefix/eigenvectors"]
            size(first_vectors, 2) == length(first_values) || error("JLD2 eigenpair shape mismatch")
            all(isfinite, first_values) && all(isfinite, first_vectors) || error("JLD2 contains nonfinite eigenpairs")
        end
    end
    true
end

function case_row(TX, scenario, portfolio, g, screen, result, spectrum_path, model_hash, settings, pf_status, pf, failure, run_id)
    rowid = string(scenario[:scenario_id])
    numeric_ok = result !== nothing
    pfgate = pf_status == "PASS"
    feasible = isempty(screen.failures)
    fcode = String[]
    !isempty(failure) && push!(fcode, failure)
    append!(fcode, screen.failures)
    if numeric_ok
        result.dae_f_residual > settings.dae_f_tolerance && push!(fcode, "DAE_F_RESIDUAL_EXCEEDED")
        result.dae_g_residual > settings.dae_g_tolerance && push!(fcode, "DAE_G_RESIDUAL_EXCEEDED")
        result.eigenpair_residual_max > settings.eigenpair_tolerance && push!(fcode, "EIGENPAIR_RESIDUAL_EXCEEDED")
    end
    unique!(fcode)
    status = !pfgate ? "FAILED" : (!numeric_ok ? "FAILED" : (!feasible ? "PHYSICALLY_INFEASIBLE" : (!isempty(fcode) ? "NUMERICAL_GATE_FAILED" : "PASS")))
    equilibrium_hash = numeric_ok ? base_case_hash(TX, result, rowid, portfolio.id, g) : ""
    (scenario_id=rowid,
     scenario_index_zero_based=parse(Int, string(scenario[:scenario_index_zero_based])),
     portfolio_id=portfolio.id, portfolio_mask=portfolio.mask, cardinality=length(portfolio.buses),
     g=g, treatment=(g == 0.03625 ? "ORIGINAL" : "INTERVENTION"), status=status,
     load_scale=parse(Float64, string(scenario[:load_scale_L_s])),
     epsilon_vector_json=string(scenario[:epsilon_final_by_bus_pu_json]),
     generator_schedule_json=string(scenario[:P_G_scheduled_by_bus_pu_json]),
     saturation_count=length(parse_int_array(string(scenario[:redispatch_saturated_buses_json]))),
     pf_status=pf_status, pf_residual=pf === nothing ? NaN : pf[2],
     power_flow_iterations=pf === nothing ? 0 : pf[3],
     voltage_min=screen.vmin, voltage_max=screen.vmax,
     dae_f_residual=numeric_ok ? result.dae_f_residual : NaN,
     dae_g_residual=numeric_ok ? result.dae_g_residual : NaN,
     state_count=numeric_ok ? result.state_count : 0,
     transverse_dimension=numeric_ok ? result.transverse_dimension : 0,
     alpha_perp=numeric_ok ? result.alpha_perp : NaN,
     alpha_omega=numeric_ok ? result.alpha_omega : NaN,
     critical_lambda_real=numeric_ok ? real(result.critical_lambda) : NaN,
     critical_lambda_imag=numeric_ok ? imag(result.critical_lambda) : NaN,
     critical_frequency_hz=numeric_ok ? result.critical_frequency_hz : NaN,
     critical_damping_ratio=numeric_ok ? result.critical_damping_ratio : NaN,
     critical_mode_id=numeric_ok ? result.critical_mode_id : "",
     eigenpair_residual=numeric_ok ? result.eigenpair_residual_max : NaN,
     complement_condition=numeric_ok ? result.gz_condition : NaN,
     failure_code=join(fcode, ";"),
     spectrum_file=spectrum_path, run_id=run_id,
     model_hash=model_hash, equilibrium_hash=equilibrium_hash,
     alpha_perp_complete_spectrum=numeric_ok && result.state_count > 0,
     alpha_omega_available=numeric_ok && isfinite(result.alpha_omega),
     pf_case_valid=pfgate, physical_feasible=feasible)
end

function append_case_rows(path::String, rows)
    isempty(rows) && return
    exists = isfile(path)
    CSV.write(path, DataFrame(rows); append=exists, writeheader=!exists)
end

function scenario_case_row(TX, scenario, portfolio, g, screen, pf, model_hash, settings,
                           spectrum_path, run_id, result, failure)
    pf_status = pf === nothing ? "FAILED" : "PASS"
    case_row(TX, scenario, portfolio, g, screen, result, spectrum_path, model_hash,
             settings, pf_status, pf, failure, run_id)
end

function run_scenario(TX, scenario, settings, model_hash, run_root, phase, attempt)
    scenario_id = string(scenario[:scenario_id])
    input_maps = load_scenario!(TX, scenario)
    pf = nothing
    screen = (failures=String[], vmin=NaN, vmax=NaN,
              gen_pq=Dict{Int,Tuple{Float64,Float64}}(), slack_p=NaN, slack_q=NaN)
    pf_failure = ""
    try
        t0 = time()
        v, residual, iterations = TX.solve_powerflow()
        elapsed = time() - t0
        isfinite(residual) || error("PF_NONFINITE")
        residual < settings.pf_residual_tolerance || error("PF_RESIDUAL_EXCEEDED")
        pf = (v, residual, iterations)
        screen = physical_screen(TX, v, input_maps.p_gen, settings)
        settings.pf_elapsed_limit > 0 && elapsed > settings.pf_elapsed_limit && (pf_failure = "TIMEOUT_POWER_FLOW")
    catch err
        message = sprint(showerror, err)
        pf_failure = occursin("PF_NONFINITE", message) ? "PF_NONFINITE" :
                     occursin("PF_RESIDUAL_EXCEEDED", message) ? "PF_RESIDUAL_EXCEEDED" :
                     occursin("TIMEOUT", message) ? "TIMEOUT_POWER_FLOW" :
                     occursin("stalled", message) || occursin("failed", message) ? "PF_NO_CONVERGENCE" :
                     "UNCLASSIFIED_EXCEPTION"
    end

    payloads = NamedTuple[]
    output_leaf = phase == "smoke" ? "smoke_attempt_$(attempt)" : "scenarios"
    spectra_dir = joinpath(run_root, "raw", "julia", output_leaf)
    mkpath(spectra_dir)
    spectrum_path = joinpath(spectra_dir, "$(scenario_id).jld2")
    spectrum_rel = relpath(spectrum_path, run_root)
    result_rows = NamedTuple[]
    cases = [(portfolio=p, g=0.03625) for p in portfolio_spec()]
    push!(cases, (portfolio=last(portfolio_spec()), g=0.25))

    for (case_index, case) in enumerate(cases)
        portfolio, g = case.portfolio, case.g
        result = nothing
        failure = pf_failure
        if pf !== nothing
            try
                result = evaluate_case(TX, portfolio.buses, g, pf)
                if result.dae_f_residual > settings.dae_f_tolerance
                    failure = "DAE_F_RESIDUAL_EXCEEDED"
                elseif result.dae_g_residual > settings.dae_g_tolerance
                    failure = "DAE_G_RESIDUAL_EXCEEDED"
                elseif result.eigenpair_residual_max > settings.eigenpair_tolerance
                    failure = "EIGENPAIR_RESIDUAL_EXCEEDED"
                end
            catch err
                message = sprint(showerror, err)
                failure = occursin("DAE_NONFINITE", message) ? "DAE_NONFINITE" :
                          occursin("EIGENPAIR_RESIDUAL_EXCEEDED", message) ? "EIGENPAIR_RESIDUAL_EXCEEDED" :
                          occursin("SPECTRUM_NONFINITE", message) ? "SPECTRUM_NONFINITE" :
                          occursin("EIGENSOLVER_FAILURE", message) ? "EIGENSOLVER_FAILURE" :
                          "UNCLASSIFIED_EXCEPTION"
            end
        end
        result !== nothing && push!(payloads, merge(result, (portfolio_id=portfolio.id, g=g, case_index=case_index)))
        push!(result_rows, scenario_case_row(TX, scenario, portfolio, g, screen, pf,
                    model_hash, settings, spectrum_rel, basename(run_root), result, failure))
    end
    write_scenario_spectrum(spectrum_path, scenario_id, payloads)
    phase == "smoke" && validate_scenario_spectrum(spectrum_path, scenario_id, length(payloads))
    result_rows, (scenario_id=scenario_id, pf_status=pf === nothing ? "FAILED" : "PASS",
                  case_count=length(result_rows), spectrum_case_count=length(payloads),
                  physical_failures=screen.failures, spectrum_path=spectrum_rel)
end

function run_main()
    run_root = abspath(arg_value("--run-root"; default=arg_value("--campaign-root"; default="")))
    manifest_path = abspath(arg_value("--manifest"; default=""))
    phase = arg_value("--phase"; default="smoke")
    phase in ("smoke", "full") || error("--phase must be smoke or full")
    isempty(run_root) && error("--run-root required")
    isempty(manifest_path) && error("--manifest required")
    mkpath(joinpath(run_root, "reports"))
    mkpath(joinpath(run_root, "src"))
    rows = load_scenarios(manifest_path)
    length(rows) == 1000 || error("manifest must contain exactly 1000 rows; observed $(length(rows))")
    [string(r[:scenario_id]) for r in rows] == [@sprintf("IAS26-060-S%04d", i) for i in 1:1000] || error("scenario IDs are not the frozen ordered sequence")

    repo = REPO_ROOT
    runner = TX4_SOURCE_PATH
    frozen_data = TX4_DATA_PATH
    TX = TX4
    isdefined(TX, :P4_G) || error("runtime-selectable P4_G not defined")

    attempt = arg_value("--attempt"; default="01")
    settings = (v_min=parse(Float64, arg_value("--v-min"; default="0.9")),
                v_max=parse(Float64, arg_value("--v-max"; default="1.1")),
                limit_tolerance=parse(Float64, arg_value("--limit-tolerance"; default="1e-8")),
                slack_p_min=parse(Float64, arg_value("--slack-p-min"; default="0.0")),
                slack_p_max=parse(Float64, arg_value("--slack-p-max"; default="11.0")),
                pf_residual_tolerance=parse(Float64, arg_value("--pf-residual-tolerance"; default="1e-9")),
                pf_elapsed_limit=parse(Float64, arg_value("--pf-timeout"; default="60")),
                dae_f_tolerance=parse(Float64, arg_value("--dae-f-tolerance"; default="1e-7")),
                dae_g_tolerance=parse(Float64, arg_value("--dae-g-tolerance"; default="1e-7")),
                eigenpair_tolerance=parse(Float64, arg_value("--eigenpair-tolerance"; default="1e-8")))
    indices = phase == "smoke" ? select_smoke_rows(rows) : rows
    result_path = phase == "smoke" ? joinpath(run_root, "derived", "MC_CASES_SMOKE_ATTEMPT_$(attempt).csv") : joinpath(run_root, "derived", "MC_CASES.csv")
    isfile(result_path) && error("refusing to overwrite existing immutable output: $result_path")
    model_hash = bytes2hex(sha256(codeunits(sha256_file(runner) * sha256_file(frozen_data))))
    println("IAS26_060_JULIA_START phase=$phase scenarios=$(length(indices)) cases_expected=$(17*length(indices)) run_id=$(basename(run_root)) model_hash=$model_hash")
    flush(stdout)
    for (index, scenario) in enumerate(indices)
        case_rows, meta = run_scenario(TX, scenario, settings, model_hash, run_root, phase, attempt)
        append_case_rows(result_path, case_rows)
        println("IAS26_060_SCENARIO_COMPLETE phase=$phase completed=$index/$(length(indices)) id=$(meta.scenario_id) pf=$(meta.pf_status) spectra=$(meta.spectrum_case_count)/17")
        flush(stdout)
    end
    println("IAS26_060_JULIA_COMPLETE phase=$phase result_csv=$(relpath(result_path, run_root))")
end

run_main()
