# T5: nonlinear delayed validation of the closed-form delay transport on IEEE-39 (method of steps, exact delay).
# Design A: all ten PLLs at 44 ms with gains transported to protect the least-damped ~4.9 Hz mode.
# Design B: same 44 ms with the ORIGINAL gains (comparator; linear count says 8 unstable roots).
# Read-only use of the frozen model code; writes only into research_gold/checks.
using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
BLAS.set_num_threads(1)

function transport(kp, ki, lam, h)
    a, w = real(lam), imag(lam); e = exp(a * h)
    (e * (kp * cos(w * h) + (ki + a * kp) / w * sin(w * h)),
     e * (ki * cos(w * h) - (abs2(lam) * kp + a * ki) / w * sin(w * h)))
end

function main()
    ctx = R.N.design_context(R.ROOT)
    base = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
    rho, kp0, ki0 = Float64.(base["rho"]), Float64.(base["Kp"]), Float64.(base["Ki"])
    lam = -0.8015607744782098 + 30.850471693237324im      # protected mode (T2), 4.91 Hz
    h = 0.004; tau = 0.044
    g = [transport(kp0[i], ki0[i], lam, h) for i in 1:10]
    kpA, kiA = first.(g), last.(g)
    println("limits kp ", extrema(ctx.kpmin), extrema(ctx.kpmax), " ki ", extrema(ctx.kimin), extrema(ctx.kimax))
    println("design A kp ", kpA[1], " ki ", kiA[1], " within limits: ",
            all(ctx.kpmin .<= kpA .<= ctx.kpmax) && all(ctx.kimin .<= kiA .<= ctx.kimax)); flush(stdout)
    rows = NamedTuple[]
    plan = vcat([("A_transport_44ms", kpA, kiA, c, 60.0) for c in DelayedEvents.CASES],
                [("B_original_gains_44ms", kp0, ki0, DelayedEvents.CASES[2], 20.0)])
    for (name, kp, ki, (bus, delta), horizon) in plan
        println("EVENT_START ", name, " ", bus, " ", delta); flush(stdout)
        row = try
            m = R.model(ctx, rho, kp, ki; bus, delta, dc_convention=:physical_supply)
            met, r = DelayedEvents.adaptive(m; tau=tau, horizon=horizon, dtmax=.01, tol=1e-9)
            pass = met.Fpeak_Hz <= .5 && met.Rpeak_Hz_s <= .5 && met.Vmin >= .9 && met.Vmax <= 1.1 && met.limiter_fraction >= .002
            (; design=name, bus, delta, horizon, F=met.Fpeak_Hz, R=met.Rpeak_Hz_s, Vmin=met.Vmin, Vmax=met.Vmax, slack=met.limiter_fraction, pass, error="")
        catch err
            (; design=name, bus, delta, horizon, F=NaN, R=NaN, Vmin=NaN, Vmax=NaN, slack=NaN, pass=false, error=first(sprint(showerror, err), 200))
        end
        push!(rows, row); CSV.write(joinpath(@__DIR__, "T5_nonlinear_events.csv"), DataFrame(rows))
        println("EVENT_DONE ", row); flush(stdout)
    end
end
main()
