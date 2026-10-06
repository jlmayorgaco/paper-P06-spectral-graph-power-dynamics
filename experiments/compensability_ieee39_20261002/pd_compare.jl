# Counterfactual PowerDynamics run: hold rho fixed and change only PLL gains.
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "experiments", "codesign_validation_20261001", "validate_pd.jl"))
include(joinpath(ROOT, "experiments", "nonlinear_codesign_20261001", "PDPhysicalReference.jl"))

using .PDPhysicalReference
using CSV, DataFrames, TOML, LinearAlgebra, SHA
using OrdinaryDiffEqRosenbrock, SciMLBase
using NetworkDynamics

const RHO_PATH = joinpath(ROOT, "reports", "nonlinear_codesign_20261001", "candidate_final_physical.toml")
const COMP_OUT = joinpath(ROOT, "reports", "poster", "ias2026", "compensability_ieee39_20261002")
mkpath(COMP_OUT)

function collect_response(sol, dt, window)
    times = collect(0.0:dt:61.0)
    phase = zeros(39, length(times))
    volt = zeros(39, length(times))
    previous = zeros(39)
    initial = zeros(39)
    unwrap = zeros(39)

    for (j, t) in enumerate(times)
        ss = NWState(sol, abs(t - 1.0) < 1e-12 ? t + 1e-9 : t)
        for b in 1:39
            ur = Float64(ss[VIndex(b, Symbol("busbar₊u_r"))])
            ui = Float64(ss[VIndex(b, Symbol("busbar₊u_i"))])
            raw = atan(ui, ur)
            if j == 1
                previous[b] = raw
                initial[b] = raw
                unwrap[b] = raw
            else
                unwrap[b] += mod(raw - previous[b] + pi, 2pi) - pi
                previous[b] = raw
            end
            phase[b, j] = unwrap[b] - initial[b]
            volt[b, j] = hypot(ur, ui)
        end
    end

    lag = round(Int, window / dt)
    frequency = zeros(size(phase))
    rocof = zeros(size(phase))
    for j in eachindex(times), b in 1:39
        a = phase[b, j]
        p1 = j > lag ? phase[b, j-lag] : 0.0
        p2 = j > 2lag ? phase[b, j-2lag] : 0.0
        frequency[b, j] = (a-p1)/(2pi*window)
        rocof[b, j] = (a-2p1+p2)/(2pi*window^2)
    end
    fi = argmax(abs.(frequency))
    ri = argmax(abs.(rocof))
    (; Fpeak_Hz=maximum(abs, frequency), Rpeak_Hz_s=maximum(abs, rocof),
       F_bus=fi[1], F_time_s=times[fi[2]]-1,
       R_bus=ri[1], R_time_s=times[ri[2]]-1,
       Vmin_pu=minimum(volt), Vmax_pu=maximum(volt))
end

function simulate_arm(label, rho, kp, ki, base)
    println("ARM_START ", label); flush(stdout)
    nw = PDPhysicalReference.build_architecture(base, rho, kp, ki)
    state = PDPhysicalReference.trim_state(nw, base, rho, kp, ki)
    trim = PDPhysicalReference.residual_audit(nw, state)
    trim.maximum < 1e-7 || error("PowerDynamics trim failed for $label: $(trim.maximum)")
    dt = 0.005
    window = 0.5
    rows = NamedTuple[]
    for bus in (8, 16, 29), delta in (100.0, -100.0)
        println("PD_START ", label, " bus=", bus, " delta_MW=", delta); flush(stdout)
        active = event_network(nw, bus, delta)
        started = time()
        watch = DiscreteCallback((u, t, int) -> time()-started > 180.0,
            int -> terminate!(int); save_positions=(false, false))
        elapsed = @elapsed sol = SciMLBase.solve(
            SciMLBase.ODEProblem(active, state, (0.0, 61.0)), Rodas5P();
            callback=CallbackSet(get_callbacks(active), watch),
            initializealg=SciMLBase.NoInit(), saveat=dt,
            abstol=1e-9, reltol=1e-9, maxiters=250000)
        complete = SciMLBase.successful_retcode(sol.retcode) && sol.t[end] >= 61.0-1e-8
        if !complete
            row = (; arm=label, bus, delta_MW=delta, complete=false,
                   last_time_s=sol.t[end], retcode=string(sol.retcode), runtime_s=elapsed)
        else
            met = collect_response(sol, dt, window)
            row = (; arm=label, bus, delta_MW=delta, complete=true,
                   last_time_s=sol.t[end], retcode=string(sol.retcode), runtime_s=elapsed,
                   met..., F_pass=met.Fpeak_Hz <= 0.5, R_pass=met.Rpeak_Hz_s <= 0.5,
                   V_pass=met.Vmin_pu >= 0.9 && met.Vmax_pu <= 1.1)
        end
        push!(rows, row)
        CSV.write(joinpath(COMP_OUT, "$(label)_events.csv"), DataFrame(rows))
        println("PD_DONE ", row); flush(stdout)
    end
    rows
end

function main()
    candidate = TOML.parsefile(RHO_PATH)
    rho = Float64.(candidate["rho"])
    kp_tuned = Float64.(candidate["Kp"])
    ki_tuned = Float64.(candidate["Ki"])
    kp_nominal = fill(PDPhysicalReference.NOMINAL_KP, 10)
    ki_nominal = fill(PDPhysicalReference.NOMINAL_KI, 10)
    base = PDPhysicalReference.frozen_baseline()
    tuned_rows = simulate_arm("tuned_PLL", rho, kp_tuned, ki_tuned, base)
    nominal_rows = simulate_arm("nominal_PLL_same_rho", rho, kp_nominal, ki_nominal, base)
    combined = vcat(tuned_rows, nominal_rows)
    CSV.write(joinpath(COMP_OUT, "events.csv"), DataFrame(combined))
    tuned_nw = PDPhysicalReference.build_architecture(base, rho, kp_tuned, ki_tuned)
    tuned_state = PDPhysicalReference.trim_state(tuned_nw, base, rho, kp_tuned, ki_tuned)
    trim = PDPhysicalReference.residual_audit(tuned_nw, tuned_state)

    open(joinpath(COMP_OUT, "protocol.toml"), "w") do io
        TOML.print(io, Dict(
            "candidate_sha256" => bytes2hex(SHA.sha256(read(RHO_PATH))),
            "rho_held_fixed_between_arms" => true,
            "tuned_gains_from_candidate" => true,
            "nominal_Kp" => PDPhysicalReference.NOMINAL_KP,
            "nominal_Ki" => PDPhysicalReference.NOMINAL_KI,
            "events" => [[b, d] for b in (8,16,29) for d in (100.0,-100.0)],
            "backend" => "compiled PowerDynamics DAE",
            "horizon_s_after_event" => 60.0,
            "sample_step_s" => 0.005,
            "monitor_buses" => collect(1:39),
            "trim_residual_tuned" => trim.maximum,
            "claim" => "Finite six-case counterfactual at fixed replacement fractions; no continuous robust or global optimality claim."
        ))
    end
    println("EXPERIMENT_DONE tuned_complete=", count(r -> get(r, :complete, false), tuned_rows),
        " nominal_complete=", count(r -> get(r, :complete, false), nominal_rows))
end

main()
