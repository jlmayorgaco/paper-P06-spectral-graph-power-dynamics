# ---------------------------------------------------------------------------
# sim_ieee39_event.jl  --  "truth" generator for the virtual-PMU experiment
# PowerDynamics.jl v5.0.0 IEEE-39 example (native physics: 10 Sauer-Pai
# machines, AVRTypeI, TGOV1, constant-Z loads), with ONE event between t0 and
# tf, after which the parameter is restored and the grid returns to normal.
#
# usage:
#   julia --project=<env> sim_ieee39_event.jl OUTDIR TYPE WHERE SCALE T0 TF TEND FPS
#     TYPE  = load | gen
#     WHERE = bus number (load bus for TYPE=load, generator bus 30..38 for gen)
#     SCALE = relative change  (load: Pset,Qset *= 1+SCALE ; gen: p_ref *= 1+SCALE)
# writes OUTDIR/truth_V.csv (t, Re/Im V of the 39 buses), OUTDIR/loads.csv,
#        OUTDIR/meta.json and a copy of the IEEE-39 CSV data used.
# ---------------------------------------------------------------------------
using PowerDynamics, PowerDynamics.Library, ModelingToolkitBase, NetworkDynamics
using NetworkDynamics: SII
using DataFrames, CSV, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve, SciMLBase
using Printf, JSON

OUT, TYPE, WHERE, SCALE, T0, TF, TEND, FPS = ARGS[1], ARGS[2], parse(Int, ARGS[3]),
    parse(Float64, ARGS[4]), parse(Float64, ARGS[5]), parse(Float64, ARGS[6]),
    parse(Float64, ARGS[7]), parse(Int, ARGS[8])
mkpath(OUT)
DATA = joinpath(pkgdir(PowerDynamics), "docs", "examples", "ieee39data")
for f in readdir(DATA); cp(joinpath(DATA, f), joinpath(OUT, f); force=true); end

# ---------------- build the network exactly as the PowerDynamics tutorial ----
branch_df = CSV.read(joinpath(DATA, "branch.csv"), DataFrame)
bus_df = CSV.read(joinpath(DATA, "bus.csv"), DataFrame)
load_df = CSV.read(joinpath(DATA, "load.csv"), DataFrame)
machine_df = CSV.read(joinpath(DATA, "machine.csv"), DataFrame)
avr_df = CSV.read(joinpath(DATA, "avr.csv"), DataFrame)
gov_df = CSV.read(joinpath(DATA, "gov.csv"), DataFrame)
set_Sbase!(100.0); set_fbase!(60.0)
load = ZIPLoad(; name=:ZIPLoad, Pset=nothing, Qset=nothing, KpZ=nothing, KqZ=nothing,
    KpI=nothing, KqI=nothing, KpC=nothing, KqC=nothing)
mp = (; Sn=nothing, Vn=nothing, R_s=nothing, X_ls=nothing, X_d=nothing, X_q=nothing,
    X′_d=nothing, X′_q=nothing, X″_d=nothing, X″_q=nothing, T′_d0=nothing, T′_q0=nothing,
    T″_d0=nothing, T″_q0=nothing, H=nothing, D=nothing)
unctrl = SauerPaiMachine(; τ_m_input=false, vf_input=false, name=:machine, mp...)
ctrl = CompositeInjector([SauerPaiMachine(; name=:machine, mp...),
    AVRTypeI(; name=:avr, ceiling_function=:quadratic, Ka=nothing, Ke=nothing, Kf=nothing,
        Ta=nothing, Tf=nothing, Te=nothing, Tr=nothing, vr_min=nothing, vr_max=nothing,
        E1=nothing, Se1=nothing, E2=nothing, Se2=nothing),
    TGOV1(; name=:gov, V_min=nothing, V_max=nothing, R=nothing, T1=nothing, T2=nothing,
        T3=nothing, DT=nothing, ω_ref=nothing)], name=:ctrld_gen)
tmpl = Dict("junction" => compile_bus(MTKBus(); name=:j),
    "load" => compile_bus(MTKBus(load); name=:l),
    "ctrld_machine" => compile_bus(MTKBus(ctrl); name=:cm),
    "ctrld_machine_load" => compile_bus(MTKBus(ctrl, load); name=:cml),
    "unctrld_machine_load" => compile_bus(MTKBus(unctrl, load); name=:uml))
function apply!(bus, tbl, i)
    row = tbl[findfirst(tbl.bus .== i), :]
    for c in names(tbl)
        (c == "bus" || c == "V_b") && continue
        set_default!(bus, Regex(c * "\$"), row[c])
    end
end
busses = map(eachrow(bus_df)) do row
    b = compile_bus(tmpl[row.category]; vidx=row.bus, name=Symbol("bus$(row.bus)"))
    set_default!(b, :busbar₊Vbase, row.base_kv)
    row.has_load && apply!(b, load_df, row.bus); row.has_gen && apply!(b, machine_df, row.bus)
    row.has_avr && apply!(b, avr_df, row.bus); row.has_gov && apply!(b, gov_df, row.bus)
    set_pfmodel!(b, row.bus_type == "PQ" ? pfPQ(P=row.P, Q=row.Q) :
                    row.bus_type == "PV" ? pfPV(P=row.P, V=row.V) : pfSlack(V=row.V, δ=0))
    b
end
lt = compile_line(MTKLine(PiLine_fault(; name=:piline)); name=:pl)
lines = map(eachrow(branch_df)) do row
    l = compile_line(lt; src=row.src_bus, dst=row.dst_bus)
    for c in names(branch_df)
        c ∉ ["src_bus", "dst_bus", "transformer"] && set_default!(l, Regex(c * "\$"), row[c])
    end
    l
end
nw = Network(busses, lines)
f = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], f); set_initformula!(nw[VIndex(39)], f)
s0 = initialize_from_pf!(nw; verbose=false)
u0 = copy(uflat(s0)); p0 = copy(pflat(s0))
prob = ODEProblem(nw, copy(u0), (0.0, TEND), copy(p0))
ir = [SII.variable_index(prob, VIndex(b, :busbar₊u_r)) for b in 1:39]
ii = [SII.variable_index(prob, VIndex(b, :busbar₊u_i)) for b in 1:39]

# ---------------- event parameters ----------------
pidx = TYPE == "load" ?
    [SII.parameter_index(prob, VIndex(WHERE, :ZIPLoad₊Pset)), SII.parameter_index(prob, VIndex(WHERE, :ZIPLoad₊Qset))] :
    [SII.parameter_index(prob, VIndex(WHERE, :ctrld_gen₊gov₊p_ref))]
params(active) = (p = copy(p0); active && (p[pidx] .*= (1 + SCALE)); p)

# ---------------- three-segment integration: [0,t0) [t0,tf) [tf,TEND] -------
tgrid = collect(0:1/FPS:TEND)
T, U = Float64[], Vector{Vector{Float64}}()
u = copy(u0)
for (k, (a, b, act)) in enumerate([(0.0, T0, false), (T0, TF, true), (TF, TEND, false)])
    ts = filter(t -> t >= a - 1e-12 && (k == 3 ? t <= b + 1e-12 : t < b - 1e-12), tgrid)
    pr = remake(prob; u0=copy(u), p=params(act), tspan=(a, b))
    sol = solve(pr, Rodas5P(); reltol=1e-8, abstol=1e-9, dtmax=1 / 120, saveat=ts,
        initializealg=BrownFullBasicInit(abstol=1e-12))   # differential states kept, algebraic re-solved
    SciMLBase.successful_retcode(sol) || error("simulation failed in segment $k: $(sol.retcode)")
    for (t, x) in zip(sol.t, sol.u)
        any(abs.(ts .- t) .< 1e-9) && (isempty(T) || t > T[end] + 1e-9) && (push!(T, t); push!(U, copy(x)))
    end
    global u = copy(sol.u[end])
end
open(joinpath(OUT, "truth_V.csv"), "w") do io
    println(io, join(["t"; ["Vr$b" for b in 1:39]; ["Vi$b" for b in 1:39]], ","))
    for (t, x) in zip(T, U)
        println(io, join([@sprintf("%.6f", t); [@sprintf("%.12g", x[i]) for i in ir]; [@sprintf("%.12g", x[i]) for i in ii]], ","))
    end
end
# initialized load parameters (the estimator's nominal load model)
open(joinpath(OUT, "loads.csv"), "w") do io
    println(io, "bus,Pset,Qset,Vset")
    for b in load_df.bus
        g(s) = p0[SII.parameter_index(prob, VIndex(b, s))]
        println(io, b, ",", g(:ZIPLoad₊Pset), ",", g(:ZIPLoad₊Qset), ",", g(:ZIPLoad₊Vset))
    end
end
open(joinpath(OUT, "meta.json"), "w") do io
    JSON.print(io, Dict("type" => TYPE, "where" => WHERE, "scale" => SCALE, "t0" => T0, "tf" => TF,
        "tend" => TEND, "fps" => FPS, "n_frames" => length(T), "powerdynamics" => string(pkgversion(PowerDynamics))))
end
println("OK frames=", length(T))
