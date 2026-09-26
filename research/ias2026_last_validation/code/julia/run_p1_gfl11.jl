#!/usr/bin/env julia

using CSV
using DataFrames
using LinearAlgebra
using Printf
using Statistics
using PowerDynamics
using NetworkDynamics

include(joinpath(@__DIR__, "frozen_gfl11.jl"))
include(joinpath(@__DIR__, "campaign_root.jl"))

const CAMPAIGN = campaign_root_from_args()
const RAW = joinpath(CAMPAIGN, "raw", "gfl11")
const REPORTS = joinpath(CAMPAIGN, "reports")
mkpath(RAW); mkpath(REPORTS)

function parse_row(row)
    x = [Float64(row[Symbol("x$(i)")]) for i in 0:10]
    f = [Float64(row[Symbol("f$(i)")]) for i in 0:10]
    return x, f
end

function relerr(a, b)
    return abs(a-b)/max(1.0, abs(a), abs(b))
end

function write_transfer(values)
    path = joinpath(RAW, "p1_transfer_julia.csv")
    open(path, "w") do io
        println(io, "frequency_hz,sigma_min,condition,h11_re,h11_im,h12_re,h12_im,h21_re,h21_im,h22_re,h22_im")
        for r in values
            println(io, join((r.frequency_hz, r.sigma_min, r.condition,
                              real(r.h11), imag(r.h11), real(r.h12), imag(r.h12),
                              real(r.h21), imag(r.h21), real(r.h22), imag(r.h22)), ','))
        end
    end
    return path
end

function main()
    oracle_path = joinpath(RAW, "p1_python_oracle.csv")
    isfile(oracle_path) || error("Run src/python/run_p1_gfl11.py before the Julia parity runner")
    df = CSV.read(oracle_path, DataFrame)
    max_state = 0.0; max_current = 0.0; max_power = 0.0
    errors = Float64[]
    for row in eachrow(df)
        x, expected_f = parse_row(row)
        f, ir, ii, ps, qs = gfl11_eval(x, row.v, row.a, row.p_ref, row.q_ref, row.v_ref, row.g, row.w)
        append!(errors, [relerr(f[i], expected_f[i]) for i in eachindex(f)])
        max_state = max(max_state, maximum(relerr.(f, expected_f)))
        max_current = max(max_current, relerr(ir, row.i_r), relerr(ii, row.i_i))
        max_power = max(max_power, relerr(ps, row.P_system), relerr(qs, row.Q_system))
    end

    # Structural PowerDynamics check: the exact transcription is a valid
    # injector, and compile_bus can construct its current-source interface.
    set_ωbase!(GFL11_OMEGA_B)
    model = GFL11Injector(name=:gfl11)
    injector_ok = isinjectormodel(model)
    compile_ok = false
    compile_message = "not attempted"
    try
        bus = compile_bus(MTKBus(model); current_source=true)
        compile_ok = !isnothing(bus)
        compile_message = "compile_bus(MTKBus(GFL11Injector); current_source=true)"
    catch err
        compile_message = sprint(showerror, err)
    end

    # One-device/infinite-bus harness. This is a device/interface gate; the
    # IEEE-39 portfolio census is reserved for P2.
    network_ok = false
    network_residual = Inf
    network_state_count = 0
    network_repeat_delta = Inf
    network_spectrum_delta = Inf
    network_spectrum_count = 0
    network_message = "not attempted"
    network_rows = NamedTuple[]
    try
        p0, q0, w0 = 0.35, -0.08, 0.75
        function build_network(shunt_conductance, vref)
            gfl_network = GFL11Injector(name=:gfl11_network, g=0.03625, w=w0,
                                         p_ref=p0/w0, q_ref=q0/w0, v_ref=vref)
            @named gfl_bus = compile_bus(MTKBus(gfl_network); current_source=true)
            set_pfmodel!(gfl_bus, pfPQ(P=p0, Q=q0; current_source=true))
            @named shunt = DynamicParallelRCShunt(R=1/shunt_conductance, B=1e-5)
            @named network_bus = compile_bus(MTKBus(shunt))
            set_pfmodel!(network_bus, pfShunt(G=shunt_conductance, B=1e-5))
            loopback = LoopbackConnection(; src=:gfl_bus, dst=:network_bus,
                                           potential=[:u_r, :u_i], flow=[:i_r, :i_i])
            @named symbolic_slack = Library.VδConstraint(V=1.0, δ=0.0)
            @named slack_bus = compile_bus(MTKBus(symbolic_slack); pf=pfSlack(V=1.0))
            @named branch = DynamicSeriesRLBranch(R=0.01, X=0.3)
            line = compile_line(MTKLine(branch); name=:gfl_to_slack,
                                src=:network_bus, dst=:slack_bus)
            @named branch_pf = PiLine(R=0.01, X=0.3)
            line_pf = compile_line(MTKLine(branch_pf); name=:gfl_to_slack_pf)
            set_pfmodel!(line, line_pf)
            return Network([gfl_bus, network_bus, slack_bus], [loopback, line])
        end
        for shunt_conductance in (0.01, 0.05, 0.10)
            # Documented PowerDynamics current-source topology: device terminal
            # -> loopback -> dynamic shunt/network bus -> dynamic PiLine ->
            # VδConstraint slack. This is intentionally not a direct loopback
            # to the slack bus.
            pf_net = build_network(shunt_conductance, 1.0)
            pf_model = powerflow_model(pf_net)
            pf_state = solve_powerflow(pf_net; pfnw=pf_model, verbose=false)
            pf_values = interface_values(pf_state)
            network_ur = VIndex(2, :busbar₊u_r)
            network_ui = VIndex(2, :busbar₊u_i)
            haskey(pf_values, network_ur) || error("network_bus u_r interface value not found: $(collect(keys(pf_values)))")
            haskey(pf_values, network_ui) || error("network_bus u_i interface value not found: $(collect(keys(pf_values)))")
            vref = hypot(pf_values[network_ur], pf_values[network_ui])
            net = build_network(shunt_conductance, vref)
            s0 = initialize_from_pf(net; verbose=false, subverbose=false, check=:none,
                                    tol=1e-6, nwtol=1e-6)
            s1 = initialize_from_pf(net; verbose=false, subverbose=false, check=:none,
                                    tol=1e-6, nwtol=1e-6)
            nstates = length(uflat(s0))
            du = zeros(Float64, nstates)
            net(du, uflat(s0), pflat(s0), 0.0)
            residual = maximum(abs, du)
            repeat_delta = isempty(uflat(s0)) ? Inf : maximum(abs.(uflat(s0) .- uflat(s1)))
            spectrum0 = sort(collect(jacobian_eigenvals(s0)); by=z -> (real(z), imag(z)))
            spectrum1 = sort(collect(jacobian_eigenvals(s1)); by=z -> (real(z), imag(z)))
            spectrum_delta = length(spectrum0) == length(spectrum1) && !isempty(spectrum0) ?
                maximum(abs.(spectrum0 .- spectrum1)) : Inf
            push!(network_rows, (shunt_conductance=shunt_conductance,
                                 state_count=nstates, residual=residual,
                                 repeat_delta=repeat_delta,
                                 spectrum_count=length(spectrum0),
                                 spectrum_delta=spectrum_delta,
                                 status=nstates >= 11 && residual <= 1e-6 &&
                                        repeat_delta <= 1e-8 && spectrum_delta <= 1e-6 ? "PASS" : "FAIL"))
            network_state_count = max(network_state_count, nstates)
            network_residual = min(network_residual, residual)
            network_repeat_delta = max(network_repeat_delta == Inf ? 0.0 : network_repeat_delta, repeat_delta)
            network_spectrum_delta = max(network_spectrum_delta == Inf ? 0.0 : network_spectrum_delta, spectrum_delta)
            network_spectrum_count = max(network_spectrum_count, length(spectrum0))
        end
        network_ok = all(r.status == "PASS" for r in network_rows)
        network_message = "GFL current-source -> LoopbackConnection -> dynamic shunt/network bus -> PiLine -> VδConstraint slack"
    catch err
        network_message = sprint(showerror, err)
    end

    sensitivity_path = joinpath(RAW, "p1_network_sensitivity.csv")
    open(sensitivity_path, "w") do io
        println(io, "shunt_conductance,state_count,residual,repeat_delta,spectrum_count,spectrum_delta,status")
        for r in network_rows
            println(io, join((r.shunt_conductance, r.state_count, r.residual,
                              r.repeat_delta, r.spectrum_count, r.spectrum_delta, r.status), ','))
        end
    end

    xeq = gfl11_initialize(1.0, 0.0, 0.35, -0.08, 0.75)
    feq, _, _, _, _ = gfl11_eval(xeq, 1.0, 0.0, xeq[3], xeq[4], 1.0, 0.03625, 0.75)
    equilibrium_residual = maximum(abs, feq)
    freqs = sort(unique(vcat(10 .^ range(-2, 2, length=401),
                             10 .^ range(log10(0.2), log10(2.0), length=201))))
    _, _, _, _, transfer = gfl11_transfer(xeq, 1.0, 0.0, xeq[3], xeq[4], 1.0, 0.03625, 0.75, freqs)
    transfer_path = write_transfer(transfer)
    cond_max = maximum(r.condition for r in transfer)
    sigma_min = minimum(r.sigma_min for r in transfer)

    parity_pass = max_state < 1e-6 && max_current < 1e-6 && max_power < 1e-6
    gate_pass = parity_pass && injector_ok && compile_ok && network_ok && equilibrium_residual < 1e-8
    report_path = joinpath(REPORTS, "P1_GFL_PARITY.md")
    open(report_path, "w") do io
        println(io, "# P1 — frozen GFL11 Julia/PowerDynamics parity")
        println(io)
        println(io, "status: ", gate_pass ? "PASS" : "STOPPED_BY_GATE")
        println(io, "evidence_class: CANONICAL_PYTHON_VS_JULIA_DEVICE_PARITY")
        println(io, "canonical_source_manifest: raw/gfl11/canonical_source_manifest.json")
        println(io, "equation_source: canonical ieee39_devices.py + 20260911_GFL_REPRODUCTION_SPEC.md section 4.2")
        println(io, "cases: ", nrow(df))
        println(io, "max_state_relative_error: ", max_state)
        println(io, "max_terminal_current_relative_error: ", max_current)
        println(io, "max_system_base_power_relative_error: ", max_power)
        println(io, "target_median_relative_error: ", median(errors))
        println(io, "target_max_relative_error: ", maximum(errors))
        println(io, "injector_interface: ", injector_ok ? "PASS" : "FAIL")
        println(io, "compile_bus_current_source: ", compile_ok ? "PASS" : "FAIL")
        println(io, "compile_message: ", compile_message)
        println(io, "infinite_bus_harness: ", network_ok ? "PASS" : "FAIL")
        println(io, "infinite_bus_residual: ", network_residual)
        println(io, "infinite_bus_state_count: ", network_state_count)
        println(io, "infinite_bus_repeat_delta: ", network_repeat_delta)
        println(io, "infinite_bus_spectrum_count: ", network_spectrum_count)
        println(io, "infinite_bus_spectrum_delta: ", network_spectrum_delta)
        println(io, "infinite_bus_message: ", network_message)
        println(io, "network_sensitivity_csv: raw/gfl11/p1_network_sensitivity.csv")
        println(io, "network_residual_preregistered_limit: 1e-6")
        println(io, "network_state_gate: >=11")
        println(io, "equilibrium_residual: ", equilibrium_residual)
        println(io, "transfer_frequency_range_hz: 0.01–100")
        println(io, "transfer_points: ", length(freqs))
        println(io, "transfer_sigma_minimum: ", sigma_min)
        println(io, "transfer_condition_maximum: ", cond_max)
        println(io, "transfer_csv: raw/gfl11/p1_transfer_julia.csv")
        println(io, "canonical_transfer_csv: raw/gfl11/p1_transfer_canonical_python.csv")
        println(io, "transfer_comparison: raw/gfl11/p1_transfer_comparison.json")
        println(io)
        println(io, "The transfer is reported over the full requested frequency range. " *
                     "The conditioning diagnostics are evidence, not a pass criterion; " *
                     "frequencies near singularity must be excluded from any claim of " *
                     "relative transfer accuracy.")
    end
    open(joinpath(RAW, "p1_julia_environment.txt"), "w") do io
        println(io, "julia=", VERSION)
        println(io, "project=env/julia/Project.toml")
        println(io, "powerdynamics=", pkgversion(PowerDynamics))
        println(io, "networkdynamics=", pkgversion(NetworkDynamics))
        println(io, "source=frozen_gfl11.jl")
    end
    println(gate_pass ? "P1_GFL_PASS" : "P1_GFL_STOPPED",
            " cases=", nrow(df), " max_state=", @sprintf("%.3e", max_state),
            " compile=", compile_ok, " eq_resid=", @sprintf("%.3e", equilibrium_residual))
    return gate_pass ? 0 : 1
end

exit(main())
