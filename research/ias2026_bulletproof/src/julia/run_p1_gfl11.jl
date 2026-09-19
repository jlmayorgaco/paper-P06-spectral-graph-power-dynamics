#!/usr/bin/env julia

using CSV
using DataFrames
using LinearAlgebra
using Printf
using Statistics
using PowerDynamics
using NetworkDynamics

include(joinpath(@__DIR__, "frozen_gfl11.jl"))

const CAMPAIGN = normpath(joinpath(@__DIR__, "..", ".."))
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
    path = joinpath(RAW, "p1_terminal_transfer.csv")
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
    network_message = "not attempted"
    try
        p0, q0, w0 = 0.35, -0.08, 0.75
        gfl_network = GFL11Injector(name=:gfl11_network, g=0.03625, w=w0,
                                     p_ref=p0/w0, q_ref=q0/w0, v_ref=1.0)
        @named gfl_bus = compile_bus(MTKBus(gfl_network); current_source=true)
        set_pfmodel!(gfl_bus, pfPQ(P=p0, Q=q0; current_source=true))
        @named slack_bus = compile_bus(MTKBus(); pf=pfSlack(V=1.0))
        # The slack bus is the infinite bus. The loopback connection makes the
        # current-source terminal voltage exactly the slack voltage.
        loopback = LoopbackConnection(; src=:gfl_bus, dst=:slack_bus,
                                       potential=[:u_r, :u_i], flow=[:i_r, :i_i])
        net = Network([gfl_bus, slack_bus], [loopback])
        s0 = initialize_from_pf(net; verbose=false, subverbose=false, check=:none,
                                tol=1e-6, nwtol=1e-6)
        network_state_count = length(uflat(s0))
        du = zeros(Float64, network_state_count)
        if !isempty(du)
            net(du, uflat(s0), pflat(s0), 0.0)
            network_residual = maximum(abs, du)
        else
            network_residual = Inf
        end
        network_ok = network_state_count >= 11 && isfinite(network_residual) && network_residual < 1e-6
        network_message = "one GFL current-source bus + exact infinite slack bus"
    catch err
        network_message = sprint(showerror, err)
    end

    xeq = gfl11_initialize(1.0, 0.0, 0.35, -0.08, 0.75)
    feq, _, _, _, _ = gfl11_eval(xeq, 1.0, 0.0, xeq[3], xeq[4], 1.0, 0.03625, 0.75)
    equilibrium_residual = maximum(abs, feq)
    freqs = 10 .^ range(-2, 2, length=81)
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
        println(io, "evidence_class: FRESH_DEVICE_LEVEL_REPRODUCTION")
        println(io, "equation_source: 20260911_GFL_REPRODUCTION_SPEC.md section 4.2")
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
        println(io, "infinite_bus_message: ", network_message)
        println(io, "equilibrium_residual: ", equilibrium_residual)
        println(io, "transfer_frequency_range_hz: 0.01–100")
        println(io, "transfer_points: ", length(freqs))
        println(io, "transfer_sigma_minimum: ", sigma_min)
        println(io, "transfer_condition_maximum: ", cond_max)
        println(io, "transfer_csv: raw/gfl11/p1_terminal_transfer.csv")
        println(io)
        println(io, "The transfer is reported over the full requested frequency range. " *
                     "The conditioning diagnostics are evidence, not a pass criterion; " *
                     "frequencies near singularity must be excluded from any claim of " *
                     "relative transfer accuracy.")
    end
    open(joinpath(RAW, "p1_julia_environment.txt"), "w") do io
        println(io, "julia=", VERSION)
        println(io, "project=", Base.active_project())
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
