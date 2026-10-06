# E6 runner: nonlinear method-of-steps delayed events with DISTRIBUTED PLL gain matrices.
# usage: julia --project=. nhop_events.jl CASE.toml LABEL [event indices]
# CASE.toml: rho(10), tau (scalar|10), Kp(10x10), KI(10x10); row i = receiving PLL, column j = measured detector.
using LinearAlgebra, CSV, DataFrames, TOML, SciMLBase, OrdinaryDiffEqRosenbrock
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const DE = DelayedEvents; const R = DE.R
BLAS.set_num_threads(1)

mat(x) = x isa AbstractMatrix ? Float64.(x) : Float64.(reduce(hcat, x)')   # TOML array-of-arrays -> row-major matrix

# Off-diagonal injection: Boff[:,j] = sum_{i!=j} e_omega_i Kp[i,j]/pll_tau_i + e_xi_i KI[i,j]
function offinjection(m, Kp, KI)
    B = zeros(length(m.x0), 10)
    for j = 1:10, i = 1:10
        i == j && continue
        B[m.gfidx[i][4], j] += Kp[i,j] / m.gp[i].pll_tau
        B[m.gfidx[i][5], j] += KI[i,j]
    end
    B
end

# Same as DE.adaptive (tau per detector allowed), plus Boff*ed and an optional trajectory save.
function adaptive_nhop(m, Boff; tau, horizon=60., dtmax=.01, tol=1e-9, savepath=nothing)
    tauv = tau isa Number ? fill(Float64(tau), 10) : Float64.(tau)
    Bd = DE.injection(m); n = length(m.x0); sols = Any[]; ends = Float64[]; xstar = copy(m.x0)
    tmin = minimum(tauv); step = tmin
    function history(t)
        t <= 0 && return xstar
        j = clamp(searchsortedfirst(ends, t-1e-12), 1, length(sols))
        t <= last(ends)+1e-9 || error("future history")
        sols[j](clamp(t, sols[j].t[1], sols[j].t[end]))
    end
    function fun!(dx, x, p, t)
        R.rhs!(dx, x, m, t)
        ed = zeros(typeof(t), 10)
        for tj in unique(tauv)                      # one detector evaluation per distinct delay
            t - tj < -1e-12 && continue
            e = DE.detector(history(max(0, t - tj)), m)
            for j = 1:10; tauv[j] == tj && (ed[j] = e[j]); end
        end
        dx .+= Bd*(ed - DE.detector(x, m)) + Boff*ed
    end
    function jac!(J, x, p, t)
        de = R.derivatives(x, m); J .= de.Fx
        _, C, _ = DE.detector(x, m, de); J .-= Bd*C
    end
    fn = ODEFunction(fun!; jac=jac!); t0 = 0.; x = copy(xstar); started = time()
    while t0 < horizon-1e-10
        t1 = min(horizon, t0 + step)
        sol = solve(ODEProblem(fn, x, (t0, t1)), Rodas5P(); reltol=tol, abstol=tol, dtmax, dense=true, save_everystep=true, maxiters=100000)
        SciMLBase.successful_retcode(sol) || error("DDE integration $(sol.retcode)")
        push!(sols, sol); push!(ends, t1); x = copy(sol.u[end]); t0 = t1
        time()-started < 400 || error("event exceeded wall limit")
    end
    ts = collect(0:.01:horizon); states = [history(t) for t in ts]
    r = (; sol=(; t=ts, u=states), ok=true, elapsed=time()-started, retcode="METHOD_OF_STEPS_RODAS5P", last_time=horizon)
    met = savepath === nothing ? R.metrics(m, r; dt=.01, window=.5, monitor_buses=collect(1:39)) :
          R.metrics(m, r; dt=.01, window=.5, monitor_buses=collect(1:39), savepath)
    met, r
end

function main()
    d = TOML.parsefile(ARGS[1]); label = ARGS[2]
    rho = Float64.(d["rho"]); Kp = mat(d["Kp"]); KI = mat(d["KI"]); tau = d["tau"]
    tau = tau isa Number ? Float64(tau) : Float64.(tau)
    @assert size(Kp) == (10,10) && size(KI) == (10,10)
    kp = diag(Kp); ki = diag(KI)
    idx = length(ARGS) >= 3 ? parse.(Int, ARGS[3:end]) : collect(1:length(DE.CASES))
    ctx = R.N.design_context(R.ROOT); rows = NamedTuple[]
    dir = joinpath(@__DIR__, "..", "raw", "events"); mkpath(joinpath(dir, "traj"))
    taums = tau isa Number ? 1000tau : 1000maximum(tau)
    for k in idx
        bus, delta = DE.CASES[k]
        println("EVENT_START ", label, " ", k, " ", bus, " ", delta); flush(stdout)
        m = R.model(ctx, rho, kp, ki; bus, delta, dc_convention=:physical_supply)
        try
            met, _ = adaptive_nhop(m, offinjection(m, Kp, KI); tau,
                savepath=joinpath(dir, "traj", "$(label)_bus$(bus)_$(Int(delta))_trajectory.csv"))
            push!(rows, (; label, tau_ms=taums, bus, delta, F=met.Fpeak_Hz, R=met.Rpeak_Hz_s, Vmin=met.Vmin, Vmax=met.Vmax, slack=met.limiter_fraction, runtime=met.runtime_s,
                pass=met.Fpeak_Hz <= .5 && met.Rpeak_Hz_s <= .5 && met.Vmin >= .9 && met.Vmax <= 1.1 && met.limiter_fraction >= .002, error=""))
        catch err
            push!(rows, (; label, tau_ms=taums, bus, delta, F=NaN, R=NaN, Vmin=NaN, Vmax=NaN, slack=NaN, runtime=NaN, pass=false, error=sprint(showerror, err)[1:min(end,100)]))
        end
        CSV.write(joinpath(dir, "events_$(label).csv"), DataFrame(rows)); println("EVENT_DONE ", last(rows)); flush(stdout)
    end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
