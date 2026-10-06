using LinearAlgebra, CSV, DataFrames, TOML, SHA, Pkg
include(joinpath(@__DIR__, "..", "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
const OUT = @__DIR__
BLAS.set_num_threads(1)
function main()
    ctx=R.N.design_context(R.ROOT)
    design=joinpath(OUT,"..","all_pll_gain_map_validation_20261003","designs","fixed.toml")
    d=TOML.parsefile(design)
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    all(ctx.kpmin .<=kp.<=ctx.kpmax) && all(ctx.kimin .<=ki.<=ctx.kimax) || error("GAIN_LIMIT")
    m0=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
    f=zeros(length(m0.x0));R.rhs!(f,m0.x0,m0,0.)
    CSV.write(joinpath(OUT,"EQUILIBRIUM.csv"),DataFrame([(;residual=norm(f,Inf),
        GFL_MW=dot(ctx.power,rho),SG_MW=dot(ctx.power,1 .-rho),
        design_sha256=bytes2hex(sha256(read(design))),julia=string(VERSION))]))
    norm(f,Inf)<=1e-7 || error("EQUILIBRIUM_PARITY")
    rows=NamedTuple[];mkpath(joinpath(OUT,"nonlinear"))
    for (bus,delta) in DelayedEvents.CASES
        println("EVENT_START fixed ",bus," ",delta);flush(stdout)
        m=R.model(ctx,rho,kp,ki;bus,delta,dc_convention=:physical_supply)
        row=try
            met,r=DelayedEvents.adaptive(m;tau=.04,horizon=60.,dtmax=.01,tol=1e-9)
            pass=met.Fpeak_Hz<=.5 && met.Rpeak_Hz_s<=.5 && met.Vmin>=.9 &&
                 met.Vmax<=1.1 && met.limiter_fraction>=.002
            R.metrics(m,r;dt=.01,window=.5,monitor_buses=collect(1:39),
                savepath=joinpath(OUT,"nonlinear","bus$(bus)_$(Int(delta))_trajectory.csv"))
            (;design="fixed",bus,delta,F=met.Fpeak_Hz,R=met.Rpeak_Hz_s,Vmin=met.Vmin,
              Vmax=met.Vmax,slack=met.limiter_fraction,runtime=met.runtime_s,pass,error="")
        catch err
            (;design="fixed",bus,delta,F=NaN,R=NaN,Vmin=NaN,Vmax=NaN,slack=NaN,
              runtime=NaN,pass=false,error=sprint(showerror,err))
        end
        push!(rows,row);CSV.write(joinpath(OUT,"TABLE_01_FIXED_GAIN_EVENTS.csv"),DataFrame(rows))
        println("EVENT_DONE ",row);flush(stdout)
    end
    println("COMPLETE all_events_pass=",all(x.pass for x in rows));flush(stdout)
end
main()
