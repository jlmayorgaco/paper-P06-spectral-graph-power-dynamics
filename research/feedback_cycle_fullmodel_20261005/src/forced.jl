# Phase K: forced response with a small sinusoidal load-admittance modulation (peak AMP_MW at bus FORCE_BUS), exact method-of-steps delay.
# usage: julia --project=. forced.jl CASE.toml FREQ_HZ LABEL [HORIZON] [AMP_MW] [BUS]   (CASE.toml: rho, Kp, Ki, tau  -- tau = 0 means no PLL delay loop, i.e. all-SG)
using LinearAlgebra, CSV, DataFrames, TOML, SciMLBase, OrdinaryDiffEqRosenbrock
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const DE = DelayedEvents
const R = DE.R
BLAS.set_num_threads(1)
case, freq, label = ARGS[1], parse(Float64, ARGS[2]), ARGS[3]
horizon = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 30.0
amp = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 5.0
bus = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 29
d = TOML.parsefile(case); rho, kp, ki = Float64.(d["rho"]), Float64.(d["Kp"]), Float64.(d["Ki"]); tau = Float64(get(d, "tau", 0.04))
ctx = R.N.design_context(R.ROOT)
m = R.model(ctx, rho, kp, ki; bus=bus, delta=0.0, dc_convention=:physical_supply)
m1 = R.model(ctx, rho, kp, ki; bus=bus, delta=amp, dc_convention=:physical_supply)
const RED0 = copy(m.net.reduced); const LIFT0 = copy(m.net.lift)
const DRED = m1.net.reduced - RED0; const DLIFT = m1.net.lift - LIFT0
setnet!(t) = (a = sin(2pi * freq * t); m.net.reduced .= RED0 .+ a .* DRED; m.net.lift .= LIFT0 .+ a .* DLIFT; nothing)
npll = any(!isempty(g) for g in m.gfidx)
B = npll ? DE.injection(m) : zeros(length(m.x0), 10)
xstar = copy(m.x0); sols = Any[]; ends = Float64[]
function history(t)
    t <= 0 && return xstar
    j = clamp(searchsortedfirst(ends, t - 1e-12), 1, length(sols))
    sols[j](clamp(t, sols[j].t[1], sols[j].t[end]))
end
function det_at(x, t)
    setnet!(t); DE.detector(x, m)
end
function fun!(dx, x, p, t)
    setnet!(t); R.rhs!(dx, x, m, t)
    if npll && tau > 0
        ed = t - tau < -1e-12 ? zeros(10) : det_at(history(max(0, t - tau)), t - tau)
        setnet!(t)
        dx .+= B * (ed - DE.detector(x, m))
    end
end
function jac!(J, x, p, t)
    setnet!(t); de = R.derivatives(x, m); J .= de.Fx
    if npll && tau > 0; _, C, _ = DE.detector(x, m, de); J .-= B * C; end
end
function tgrad!(dT, x, p, t)
    h = 1e-6; f1 = similar(x); f2 = similar(x); fun!(f1, x, p, t + h); fun!(f2, x, p, t - h); dT .= (f1 .- f2) ./ (2h); setnet!(t); nothing
end
fn = ODEFunction(fun!; jac=jac!, tgrad=tgrad!); t0 = 0.0; x = copy(xstar); status = "ok"; started = time()
try
    global t0, x
    while t0 < horizon - 1e-10
        t1 = min(horizon, tau > 0 ? t0 + tau : horizon)
        sol = solve(ODEProblem(fn, x, (t0, t1)), Rodas5P(); reltol=1e-9, abstol=1e-9, dtmax=0.01, dense=true, save_everystep=true, maxiters=200000)
        SciMLBase.successful_retcode(sol) || error("retcode $(sol.retcode)")
        push!(sols, sol); push!(ends, t1); x = copy(sol.u[end]); t0 = t1
        time() - started < 1500 || error("wall limit")
    end
catch err
    global status = "aborted: " * sprint(showerror, err)[1:min(end, 120)]
end
tend = isempty(ends) ? 0.0 : last(ends)
ts = collect(0:0.01:tend); n = length(ts)
ebuf = zeros(n); wbuf = zeros(n); fbuf = zeros(n); vmin = Inf; vmax = -Inf
phase = zeros(10, n); prev = zeros(10); uw = zeros(10); init = zeros(10)
slackmin = Inf
for (j, t) in enumerate(ts)
    xx = history(t); setnet!(t)
    v = R.voltage(xx, m; allbus=true); vm = hypot.(v[1:2:end], v[2:2:end]); global vmin = min(vmin, minimum(vm)); global vmax = max(vmax, maximum(vm))
    raw = atan.(v[2*30:2:2*39], v[2*30-1:2:2*39-1])
    if j == 1; init .= raw; prev .= raw; end
    uw .+= mod.(raw - prev .+ pi, 2pi) .- pi; prev .= raw; phase[:, j] .= uw
    if npll
        e = DE.detector(xx, m); ebuf[j] = maximum(abs.(e) ./ vm[30:39])
        wbuf[j] = maximum(abs(xx[m.gfidx[i][4]] - xstar[m.gfidx[i][4]]) for i in 1:10 if !isempty(m.gfidx[i]))
    end
    lag = round(Int, 0.5 / 0.01)
    fbuf[j] = j > lag ? maximum(abs.(phase[:, j] - phase[:, j-lag])) / (2pi * 0.5) : 0.0
    for i in 1:10
        p = m.sp[i]; (!isempty(m.sgidx[i]) && p.controlled) || continue; ix = m.sgidx[i]
        for (jj, lo, hi) in ((1, p.gov_vmin, p.gov_vmax), (5, p.avr_vr_min, p.avr_vr_max))
            w = hi - lo; global slackmin = min(slackmin, (xx[ix[jj]] - lo) / w, (hi - xx[ix[jj]]) / w)
        end
    end
end
deg = 180 / pi; thr = 5 / deg
first = findfirst(>(thr), ebuf)
out = DataFrame(label=[label], freq_hz=[freq], amp_MW=[amp], bus=[bus], tau_ms=[1000tau], status=[status], t_end_s=[tend],
    peak_pll_err_deg=[maximum(ebuf) * deg], peak_pll_omega_dev=[maximum(wbuf)], peak_busfreq_dev_Hz=[maximum(fbuf)],
    peak_pll_err_deg_first10s=[maximum(ebuf[ts .<= 10]) * deg], peak_busfreq_dev_Hz_first10s=[maximum(fbuf[ts .<= 10])],
    t_pll_err_gt_5deg_s=[first === nothing ? NaN : ts[first]], vmin=[vmin], vmax=[vmax], sg_slack_min=[slackmin])
dir = joinpath(@__DIR__, "..", "raw", "forced"); mkpath(dir)
CSV.write(joinpath(dir, "forced_$(label).csv"), out)
sel = 1:5:n; CSV.write(joinpath(dir, "trace_$(label).csv"), DataFrame(t=ts[sel], pll_err_deg=ebuf[sel] .* deg, pll_omega_dev=wbuf[sel], busfreq_dev_Hz=fbuf[sel]))
println(out)
