# T9: parametrised nonlinear delayed events (method of steps, exact delay) for a retuning rule.
# usage: julia --project=. t9_events.jl LABEL RULE NODE_RE NODE_IM H_MS RHO HORIZON [EVENT_INDICES...]
#   RULE = exact | lowfreq | none      RHO = uniform SG->GFL share (baseline 0.875)
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
    label, rule = ARGS[1], ARGS[2]
    lam = complex(parse(Float64, ARGS[3]), parse(Float64, ARGS[4]))
    h = parse(Float64, ARGS[5]) / 1000; rhou = parse(Float64, ARGS[6]); horizon = parse(Float64, ARGS[7])
    idx = length(ARGS) > 7 ? parse.(Int, ARGS[8:end]) : collect(1:5)
    ctx = R.N.design_context(R.ROOT)
    base = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
    kp0, ki0 = Float64.(base["Kp"]), Float64.(base["Ki"]); rho = fill(rhou, 10); tau = 0.040 + h
    kp, ki = copy(kp0), copy(ki0)
    if rule == "exact"
        g = [transport(kp0[i], ki0[i], lam, h) for i in 1:10]; kp, ki = first.(g), last.(g)
    elseif rule == "lowfreq"
        kp = kp0 .+ h .* ki0
    end
    println(label, " tau=", tau, " kp=", kp[1], " ki=", ki[1], " within limits: ",
            all(ctx.kpmin .<= kp .<= ctx.kpmax) && all(ctx.kimin .<= ki .<= ctx.kimax)); flush(stdout)
    rows = NamedTuple[]
    for (bus, delta) in DelayedEvents.CASES[idx]
        row = try
            m = R.model(ctx, rho, kp, ki; bus, delta, dc_convention=:physical_supply)
            met, r = DelayedEvents.adaptive(m; tau=tau, horizon=horizon, dtmax=.01, tol=1e-9)
            pass = met.Fpeak_Hz <= .5 && met.Rpeak_Hz_s <= .5 && met.Vmin >= .9 && met.Vmax <= 1.1 && met.limiter_fraction >= .002
            (; design=label, tau_ms=1000tau, rho=rhou, kp=kp[1], ki=ki[1], bus, delta, horizon, F=met.Fpeak_Hz, R=met.Rpeak_Hz_s, Vmin=met.Vmin, Vmax=met.Vmax, slack=met.limiter_fraction, pass, error="")
        catch err
            (; design=label, tau_ms=1000tau, rho=rhou, kp=kp[1], ki=ki[1], bus, delta, horizon, F=NaN, R=NaN, Vmin=NaN, Vmax=NaN, slack=NaN, pass=false, error=first(sprint(showerror, err), 120))
        end
        push!(rows, row); CSV.write(joinpath(@__DIR__, "T9_" * label * ".csv"), DataFrame(rows))
        println("EVENT_DONE ", row); flush(stdout)
    end
end
main()
